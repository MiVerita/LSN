#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"

using namespace std;

void ReadPrimes(const string& filename, int& p1, int& p2) {
   ifstream filePrimes(filename);                      // Apre file in lettura
   if (!filePrimes) {                                  // Controlla apertura
       cerr << "Errore apertura " << filename << endl;
       exit(1);                                        // Termina in caso di errore
   }
   filePrimes >> p1 >> p2;                             // Estrae due interi primi
   filePrimes.close();                                 // Chiude il file
}
void InitializeRandom(Random& rnd, const string& seedFile, int p1, int p2) {
   ifstream fileSeed(seedFile);
   string keyword;
   int seed[4];                                        // Array per 4 semi interi
   
   // Scorre file fino a trovare la parola "RANDOMSEED"
   while(fileSeed >> keyword) {
       if(keyword == "RANDOMSEED") {
           fileSeed >> seed[0] >> seed[1] >> seed[2] >> seed[3];  // Legge i semi
           rnd.SetRandom(seed, p1, p2);                         // Imposta RNG con semi e primi
           break;                                                // Esce dal loop
       }
   }
   fileSeed.close();                                // Chiude file seme
}
vector<double> get_angle(Random& rnd, int ndata){
   vector<double>angles;
   double x, y; //theta numero campionato uniformemente tra 0 e pi
   bool flag = true;
   for(int i = 0; i < ndata; i++){
      while(flag){   
         x = rnd.Rannyu(-1,1);
         y = rnd.Rannyu();
         if(x*x + y*y < 1){
            flag = false;
            angles.push_back(acos(x/sqrt(x*x + y*y)));
         }
      }
      flag = true;
   }
   return angles;
}

int main (int argc, char *argv[]){
   Random rnd;
   int p1, p2;
   ReadPrimes("Primes", p1, p2);                     
   InitializeRandom(rnd, "seed.out", p1, p2);

   int ndata = 100000;
   int nblocks = 100;
   int data_per_block = ndata/nblocks;

   double L = 1.0;   // lunghezza dell'ago
   double d = 1.5;   // distanza tra le linee
   int hit = 0;
   vector<double>sample_angles = get_angle(rnd, ndata);

   vector<double>block_means;
   block_means.resize(nblocks);
   
   for(int i = 0; i < nblocks; i++) {
      for(int j = 0; j < data_per_block; ++j) {
         double angle = sample_angles[i*data_per_block + j];
         double x0 = rnd.Rannyu(0, d);
         if(x0 + (L/2)*sin(angle) >= d or x0 - (L/2)*sin(angle) <= 0){
            hit++;
         }
      }
      block_means[i] = 2*L*data_per_block / (hit * d);
      hit = 0;
   }

   //Calcolo delle media progressive
   double cum_mean  = 0.0;
   double cum_mean2= 0.0;
  
   ofstream results_buffon;
  
   results_buffon.open("../OUTPUT/results_buffon.dat");
  
   for(int b = 0; b < nblocks; ++b) {
      cum_mean  += block_means[b];
      cum_mean2 += block_means[b] * block_means[b];
  
      // numero di blocchi finora
      int m = b + 1;
  
      // media progressiva fino al blocco b
      double mean_b = cum_mean  / m;
      double mean2_b = cum_mean2 / m;
  
      double err_mean = (b == 0)
                        ? 0.0
                        : sqrt((mean2_b - mean_b * mean_b) / (m - 1));
  
      results_buffon << m << " " << mean_b << " " << err_mean <<endl;
          
   }
   results_buffon.close();
   rnd.SaveSeed();

   return 0;
}

