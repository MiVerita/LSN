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

   ofstream output_lattice_walk;
   ofstream output_continuum_walk;
   output_lattice_walk.open("../OUTPUT/results_lattice_walk");
   output_continuum_walk.open("../OUTPUT/results_continuum_walk");
   
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

   double block_displ_disc[100]  = {0.0}; //Stima dello spostamento per percoso su reticolo nel blocco -- Ogni cella è uno step del RW
   double block_displ_cont[100]  = {0.0}; //Stima dello spostamento per percoso nel continuo nel blocco -- Ogni cella è uno step del RW
   double best_disc[100] = {0.0}; //Media della stima di blocco dello spostamento -- Discreto
   double best_cont[100] = {0.0}; //Media della stima di blocco dello spostamento -- Continuo
   double sigma2_disc[100]  = {0.0}; //Quadrato della media blocco -- Discreto
   double sigma2_cont[100] = {0.0}; //Quadrato della media blocco -- Continuo

   double x_cont, y_cont, z_cont;
   double x_disc, y_disc, z_disc;

   int N_RW = throws_RW/blocks_RW;

   if(throws_RW % blocks_RW != 0){
      cerr << "Throws number must be a multiple of block" <<endl;
   }

   for(int i = 0; i < N_RW; i++){ //Ciclo dei blocchi
      for(int j = 0; j < steps; j++){
         block_displ_cont[j] = 0;
         block_displ_disc[j] = 0;
      }
      for(int j = 0; j < blocks_RW; j++){ //Ciclo interno ai blocchi
         x_disc = 0;
         y_disc = 0;
         z_disc = 0;

         x_cont = 0;
         y_cont = 0;
         z_cont = 0;

         for(int k = 0; k < steps; k++){ //Ciclo su random walk
            //Simulazione per cammino su reticolo
            double dim = static_cast<int>(rnd.Rannyu(1.0,4.0)); //dimensione nella quale procedere
            double dir = 2*static_cast<int>(rnd.Rannyu(0.0, 2.0)) - 1; // verso nella quale procedere

            if(dim == 1){
               x_disc += a*dir;
            }else if(dim == 2){
               y_disc += a*dir;
            }else{
               z_disc += a*dir;
            }
            double r2_disc = pow(x_disc,2) + pow(y_disc,2) + pow(z_disc,2); //Distanza quadra all'i-esimo step su reticolo
            block_displ_disc[k] += r2_disc;

            //Simulazione per cammino nel continuo
            double u = rnd.Rannyu();
            double v = rnd.Rannyu();

            double theta = acos(1.0 - 2.0 * u); //parametri sferici di incremento
            double phi = 2.0 * M_PI * v;

            x_cont += a*sin(theta)*cos(phi); //incremento spaziale
            y_cont += a*sin(theta)*sin(phi);
            z_cont += a*cos(theta);

            double r2_cont = pow(x_cont, 2) + pow(y_cont, 2) + pow(z_cont, 2); //Distanza quadra dell'i-esimo step nel continuo
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