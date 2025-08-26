#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "random.h"

using namespace std;

//this function evaluates the density function given by the psi-test -- Note -- It ignores the normalization because the metropolis alogithm erase its contribution 
double p(double x, double mu, double sigma) {
      double term1 = exp(-pow(x - mu, 2) / (2 * sigma * sigma));
      double term2 = exp(-pow(x + mu, 2) / (2 * sigma * sigma));
      return pow(term1 + term2, 2);
}
//this function evaluate the integrand given by the application of the hamiltonian operator over the psi -- Note -- the results is given in naural unit of hbar = 1, m = 1
double integrand(double x, double mu, double sigma){
   double term1 = exp(-pow(x - mu, 2) / (2 * sigma * sigma));
   double term2 = exp(-pow(x + mu, 2) / (2 * sigma * sigma));
       
   // Derivata seconda analitica di Ψ_T divisa per Ψ_T:
   double numerator = ((pow(x - mu, 2)/(pow(sigma,4)) - 1.0/(sigma*sigma)) * term1)
                    + ((pow(x + mu, 2)/(pow(sigma,4)) - 1.0/(sigma*sigma)) * term2);
   double denominator = term1 + term2;
   double kinetic = -0.5 * (numerator / denominator);
       
   // Potenziale V(x) = x^4 - (5/2) x^2
   double potential = pow(x,4) - (5.0/2.0)*pow(x,2);
       
   return kinetic + potential;
}
//this function return a vector with number sampled from the distribution p_x after an equilibration of 5000 steps
vector<double> metropolis(Random rnd, int samples, double delta, double mu, double sigma){
   vector<double> numbers;
   int equilibration = 5000;
   int total_steps = samples + equilibration;
   double acceptance = 0.;

   double x = 0.;

   for (int step = 0; step < total_steps; step++) {
      double x_new = x + delta*(2.0*rnd.Rannyu() - 1.0);
      double A = min(1.0, p(x_new,mu,sigma)/p(x,mu,sigma));
      if (rnd.Rannyu() < A) {
         x = x_new;        // passo accettato
         acceptance++;
      }
      // altrimenti x rimane invariato (rifiuto)
      if (step >= equilibration){
         numbers.push_back(x);
      }
   }
   acceptance /= total_steps;
   numbers.push_back(acceptance);
   return numbers;
}

int main(){
   Random rnd;
   int seed[4];
   int p1, p2;

   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.out");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
   
   ofstream results;
   results.open("../OUTPUT/results_integral_metropolis.dat");

   if(!results.is_open()){
      cerr << "PROBLEM: Unable to open results_results_integral_metropolis" << endl;
      return 0;
   }
   //simulation parameters
   int throws = 100000;
   int blocks = 1000;
   
   //psi-test parameters
   double sigma = 0.59;
   double mu    = 0.8;   
   double delta = 2.5;


   //support variables
   double block_energy = 0;
   double energy = 0;
   double acceptance = 0;
   double sigma_2_energy = 0;
   int N;
   int attempt = 0;

   vector<double> metro_p_x = metropolis(rnd, throws, delta, mu, sigma);
   acceptance = metro_p_x[throws];
   metro_p_x.pop_back();
   results << "ACCEPTANCE: " << acceptance<<endl;
   results << "###" <<endl;

   for(int i = 0; i < throws; i++){
      if(attempt % blocks == 0 and attempt/blocks != 0){
         energy += (block_energy/blocks); //sum of the mean in the current block
         
         N = attempt/blocks; //current block iteration
   
         sigma_2_energy += pow(block_energy/blocks, 2); //sum of the square of the blocks' mean
   
         results << energy/N << " ";
   
         if(attempt / blocks == 1){
            results << 0 << " " << N <<endl;
   
         }else{
            results << sqrt(((sigma_2_energy/N) - pow(energy/N,2) )/(N -1)) << " " << N <<endl;
         }
         block_energy = 0;
      }
      
      block_energy += integrand(metro_p_x[i], mu, sigma);
      attempt++;
   }

   results.close();
   rnd.SaveSeed();
   return 0;
}