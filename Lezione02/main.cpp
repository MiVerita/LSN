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
   output_uniform.open("OUTPUT/results_integral_uniform");   
   output_importance.open("OUTPUT/results_integral_importance");

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


   int attempt = 0; //current block iteration
   double x;

   double integral_uniform = 0.; //mean of the means up to the current block
   double integral_importance = 0.;
   double block_integral_uniform = 0.; //mean of the current block
   double block_integral_importance = 0;

   double sigma_integral_uniform = 0.; //variable to get dev std in the current block
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

      x = 1 - sqrt(1 - rnd.Rannyu()); //sampling of the density function p(x) = 2(1-x) in [0,1]
      block_integral_importance += ((M_PI/2)*cos(M_PI*x/2))/(2*(1-x));

      attempt++;
   }

   output_uniform.close();
   output_importance.close();

   //-----------------------
   ofstream output_lattice_walk;
   ofstream output_continuum_walk;
   output_lattice_walk.open("OUTPUT/results_lattice_walk");
   output_continuum_walk.open("OUTPUT/results_continuum_walk");
   
   if (!output_lattice_walk.is_open()){
      cerr << "PROBLEM: Unable to open results_lattice_walk" << endl;
      return 0;
   }
   if(!output_continuum_walk.is_open()){
      cerr << "PROBLEM: Unable to open results_continuum_walk" << endl;
      return 0;
   }

   int throws_RW  = 10000;
   int blocks_RW  = 100;
   double a = 1.0;
   int steps = 100;

   double block_displ_disc[100]  = {0.0}; //displacement estimantion of the lattice walk in block
   double block_displ_cont[100]  = {0.0}; //displacement estimation of the continuum walk in block
   double best_disc[100] = {0.0}; //mean of the blocks' displacement estimation -- discrete
   double best_cont[100] = {0.0}; //mean of the blocks' displacement estimation -- continuum
   double sigma2_disc[100]  = {0.0}; //square of the blocks' mean -- discrete
   double sigma2_cont[100] = {0.0}; //square of the blocks' mean -- continuum

   double x_cont, y_cont, z_cont;
   double x_disc, y_disc, z_disc;

   int N_RW = throws_RW/blocks_RW;

   if(N%1 != 0){
      cerr << "Throws number must be a multiple of block" <<endl;
   }

   for(int i = 0; i < N_RW; i++){ //blocks cycle
      for(int j = 0; j < steps; j++){
         block_displ_cont[j] = 0;
         block_displ_disc[j] = 0;
      }
      for(int j = 0; j < blocks_RW; j++){ //internal block cycle
         x_disc = 0;
         y_disc = 0;
         z_disc = 0;

         x_cont = 0;
         y_cont = 0;
         z_cont = 0;

         for(int k = 0; k < steps; k++){ //RW cycle
            //lattice walk simulation
            double dim = static_cast<int>(rnd.Rannyu(1.0,4.0));
            double dir = 2*static_cast<int>(rnd.Rannyu(0.0, 2.0)) - 1;

            if(dim == 1){
               x_disc += a*dir;
            }else if(dim == 2){
               y_disc += a*dir;
            }else{
               z_disc += a*dir;
            }
            double r2_disc = pow(x_disc,2) + pow(y_disc,2) + pow(z_disc,2);
            block_displ_disc[k] += r2_disc;

            //continuunm RW simulation
            double u = rnd.Rannyu();
            double v = rnd.Rannyu();

            double theta = acos(1.0 - 2.0 * u);
            double phi = 2.0 * M_PI * v;

            x_cont += a*sin(theta)*cos(phi);
            y_cont += a*sin(theta)*sin(phi);
            z_cont += a*cos(theta);

            double r2_cont = pow(x_cont, 2) + pow(y_cont, 2) + pow(z_cont, 2);
            block_displ_cont[k] += r2_cont;
         
         }
      }
      for(int j = 0; j < steps; j++){
         best_disc[j] += sqrt(block_displ_disc[j]/blocks_RW);
         sigma2_disc[j] += (block_displ_disc[j]/blocks_RW);

         best_cont[j] += sqrt(block_displ_cont[j]/blocks_RW);
         sigma2_cont[j] += (block_displ_cont[j]/blocks_RW);
      }
   }
   
   for(int i = 0; i < steps; i++){
      output_lattice_walk << best_disc[i]/N_RW << " " << sqrt(((sigma2_disc[i]/N_RW) - pow((best_disc[i]/N_RW),2))/(N_RW - 1)) << " " << i+1 <<endl;
      output_continuum_walk << best_cont[i]/N_RW << " " << sqrt(((sigma2_cont[i]/N_RW) - pow((best_cont[i]/N_RW),2))/(N_RW - 1)) << " " << i+1 <<endl;
   }

   output_lattice_walk.close();
   output_continuum_walk.close();
   
   rnd.SaveSeed();
   return 0;
}