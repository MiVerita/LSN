#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <cmath>
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

int main (int argc, char *argv[]){
   Random rnd;
   int p1, p2;
   ReadPrimes("Primes", p1, p2);                     
   InitializeRandom(rnd, "seed.out", p1, p2);
   //Definizione del numero di generazioni totali e del numero di blocchi
   int ndata = 1000000;
   int data_per_block = 10000;
   int nblocks = ndata/data_per_block;

   //Calcolo delle media di blocco
   vector<double> data(ndata);

   // generate invoca la lambda n volte pari alla dimensione del vector, assegnandone il risultato agli elementi del vettore
   generate(data.begin(), data.end(),
       [&]() {
         return rnd.Rannyu();
       }
   );

   vector<double>block_means_mean, block_means_var;
   block_means_mean.resize(nblocks);
   block_means_var.resize(nblocks);
  
   //Data blocking
   for(int i = 0; i < nblocks; i++) {
      double sum_mean = 0.0;
      double sum_var = 0.0;
      for(int j = 0; j < data_per_block; ++j) {
         sum_mean += data[i*data_per_block + j];
         sum_var += pow(data[i*data_per_block + j] - 0.5,2);
      }
      block_means_mean[i] = sum_mean / data_per_block;
      block_means_var[i] = sum_var / data_per_block;
   }

   //Calcolo delle media progressive
   double cum_mean_mean  = 0.0;
   double cum_mean2_mean = 0.0;
   double cum_mean_var  = 0.0;
   double cum_mean2_var = 0.0;
  
   ofstream results_mean,results_var;
  
   results_mean.open("../OUTPUT/results_mean.dat");
   results_var.open("../OUTPUT/results_var.dat");
  
   for(int b = 0; b < nblocks; ++b) {
      cum_mean_mean  += block_means_mean[b];
      cum_mean2_mean += block_means_mean[b] * block_means_mean[b];
      cum_mean_var  += block_means_var[b];
      cum_mean2_var += block_means_var[b] * block_means_var[b];
  
      // numero di blocchi finora
      int m = b + 1;
  
      // media progressiva fino al blocco b
      double mean_b_mean  = cum_mean_mean  / m;
      double mean2_b_mean = cum_mean2_mean / m;
      double mean_b_var  = cum_mean_var  / m;
      double mean2_b_var = cum_mean2_var / m;
  
      double err_mean = (b == 0)
                        ? 0.0
                        : sqrt((mean2_b_mean - mean_b_mean * mean_b_mean) / (m - 1));
      double err_var = (b == 0)
                        ? 0.0
                        : sqrt((mean2_b_var - mean_b_var * mean_b_var) / (m - 1));
  
      results_mean << m << " " << mean_b_mean << " " << err_mean <<endl;
      results_var << m << " " << mean_b_var << " " << err_var <<endl;
          
   }
   results_mean.close();
   results_var.close();

   
   ofstream results_chi_histo;
   results_chi_histo.open("../OUTPUT/results_chi_histo.dat"); 
   vector<int>occurencies(100); //vettore che contiene le occorrenze relative dei dati in ogni bin unitario tra 0 e 1
   vector<double>local_chi; //vettore per il salvataggio dei valori di chi quadro istantanei

   fill(occurencies.begin(), occurencies.end(), 0);
   double chi = 0; //valore istantaneo di chi quadro

   int data_per_local_chi = ndata/100;

   for(int i = 1; i <= 100; i++){ //Calcolo 100 distribuzioni e 100 valori di chi quadro istantanei
      for(int j = 0; j < data_per_local_chi; j++){
         occurencies[static_cast<int>(floor(data[data_per_local_chi*(i-1) + j]*100))]++; //uso il set di dati usato all'inizio
         //Estrae un numero, viene spostato nell'intervallo [0,100], viene presa la parte intera e viene usata per come indice del bin da incrementare
      }

      for(int k: occurencies){
         chi += pow((k - (data_per_block/nblocks)),2)/(data_per_block/nblocks); //Calcolo del chi quadro
      }
      local_chi.push_back(chi);
      results_chi_histo << chi <<endl;
      chi = 0;
      fill(occurencies.begin(), occurencies.end(), 0);
   }
   results_chi_histo.close();


   ofstream results_chi_graph;
   results_chi_graph.open("../OUTPUT/results_chi_graph.dat");
   double mean = 0;
   double mean_2 = 0;
   for(int i = 1; i <= 100; i++){
      mean+=local_chi[i-1];
      mean_2+=local_chi[i-1]*local_chi[i-1];
      double err = (i == 1) ? 0. : sqrt((mean_2/i - pow(mean/i,2)) / i-1);
      results_chi_graph << i << " " << mean/i << " " << err <<endl;
   }
   results_chi_graph.close();

   rnd.SaveSeed();

   return 0;
}

