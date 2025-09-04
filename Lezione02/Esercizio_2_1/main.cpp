#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "random.h"

using namespace std;

int main() {
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

   ofstream output_uniform; 
   ofstream output_importance; 
   output_uniform.open("../OUTPUT/results_integral_uniform");   
   output_importance.open("../OUTPUT/results_integral_importance");

   if (!output_uniform.is_open()){
      cerr << "PROBLEM: Unable to open results_uniform" << endl;
      return 0;
   }
   if(!output_importance.is_open()){
      cerr << "PROBLEM: Unable to open results_importance" << endl;
      return 0;
   }

   int throws = 100000;
   int blocks = 1000;

   if(throws % blocks != 0){
      cout << "throws must be a multiple of blocks" << endl;
      return 0;
   }


   int attempt = 0;
   double x;

   double integral_uniform = 0.; //Media dell medie di blocco corrente
   double integral_importance = 0.;
   double block_integral_uniform = 0.; //Media del blocco corrente
   double block_integral_importance = 0;

   double sigma_integral_uniform = 0.;
   double sigma_integral_importance = 0.;
   double sigma2_integral_uniform = 0.;
   double sigma2_integral_importance = 0.;

   int N;

   for(int i = 0; i < throws; i++){
      if(attempt % blocks == 0 and attempt / blocks != 0){
         N = attempt/blocks;

         integral_uniform += block_integral_uniform/blocks;
         integral_importance += block_integral_importance/blocks;

         sigma2_integral_uniform += pow(block_integral_uniform/blocks,2);
         sigma_integral_uniform = integral_uniform/N;
         sigma2_integral_importance += pow(block_integral_importance/blocks,2);
         sigma_integral_importance = integral_importance/N;
         
         output_uniform << integral_uniform/N << " ";
         output_importance << integral_importance/N << " ";

         if(N == 1){
            output_uniform << 0 << " " << N <<endl;
            output_importance << 0 << " "<< N <<endl;
         }else{
            output_uniform << sqrt(((sigma2_integral_uniform/N) - pow(sigma_integral_uniform,2) )/(N -1)) << " " << N <<endl;
            output_importance << sqrt(((sigma2_integral_importance/N) - pow(sigma_integral_importance,2) )/(N -1)) << " " << N <<endl;
         }
         block_integral_uniform = 0;
         block_integral_importance = 0;
      }
      block_integral_uniform += (M_PI/2)*cos(rnd.Rannyu()*M_PI/2);

      x = 1 - sqrt(1 - rnd.Rannyu()); //campionamento della densità di probabilità p(x) = 2(1-x) in [0,1]
      block_integral_importance += ((M_PI/2)*cos(M_PI*x/2))/(2*(1-x));

      attempt++;
   }

   output_uniform.close();
   output_importance.close();

   rnd.SaveSeed();
   return 0;
}