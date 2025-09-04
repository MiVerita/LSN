#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cstring>
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

   ofstream results_central_limit_standard, results_central_limit_exp, results_central_limit_lorentz;
   results_central_limit_standard.open("../OUTPUT/results_central_limit_standard.dat");
   results_central_limit_exp.open("../OUTPUT/results_central_limit_exp.dat");
   results_central_limit_lorentz.open("../OUTPUT/results_central_limit_lorentz.dat");

   int ndata = 100000;
   double limit_variable_standard = 0;
   double limit_variable_exponential = 0;
   double limit_variable_lorentzian = 0;

   int indexes[4] = {1,2,10,100};
   double lambda = 1.15;
   double mu = 0.0;
   double gamma = 1.0;

   //l'output su file è formattato in modo da avere un numero di righe pari a throws ed ogni riga i primi 3 numeri corrispondonp
   //ad N = 1, i successivi 3 ad N = 2 ed i successivi 6 ad N = 10 ed N = 100

   for(int i = 0; i < ndata; i++){
      for(int j: indexes){
         for(int k = 0; k < j; k++){
            limit_variable_standard += rnd.Rannyu(1,6);
            limit_variable_exponential += rnd.Exponential(lambda);
            limit_variable_lorentzian += rnd.Lorentzian(mu, gamma);
         }
         results_central_limit_standard << j << " " << limit_variable_standard/j <<endl;
         results_central_limit_exp << j << " " << limit_variable_exponential/j <<endl;
         results_central_limit_lorentz << j << " " << limit_variable_lorentzian/j <<endl;
      }

      limit_variable_standard = 0;
      limit_variable_exponential = 0;
      limit_variable_lorentzian = 0;
   }

   results_central_limit_standard.close();
   results_central_limit_exp.close();
   results_central_limit_lorentz.close();
   
   rnd.SaveSeed();
   return 0;
}

