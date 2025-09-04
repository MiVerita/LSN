#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cmath>
#include "random.h"
#include <vector>

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
//flag = true --> call
//flag = false --> put
vector<double> get_direct(bool flag, Random& rnd, double ndata, double sigma, double T, double S0, double K, double r){
    vector<double> direct_data;
    for(int i = 0; i < ndata; i++){
        double Z = rnd.Gauss(0., 1.);
        double S_T = S0*exp((r - 0.5*sigma*sigma)*T + sigma*Z*sqrt(T));
    
        double profit = (flag)
                        ? max(0., S_T - K)
                        : max(0., K - S_T);
        double discounted_profit = exp(-r*T)*profit;
        direct_data.push_back(discounted_profit);
    }

    return direct_data;
}
vector<double> get_discrete(bool flag, Random& rnd, double ndata, double sigma, double T, double nsteps, double S0, double K, double r){
    vector<double> discrete_data;
    double dT = T / nsteps;

    for(unsigned int i = 0; i < ndata; i++){
        double S_T = S0;
        for(unsigned int j = 0; j < nsteps; j++){
            double Z = rnd.Gauss(0., 1.);
            S_T *= exp((r - 0.5*sigma*sigma)*dT + sigma*Z*sqrt(dT));
        }
        double profit = (flag)
                        ? max(0., S_T - K)
                        : max(0., K - S_T);
        double discounted_profit = exp(-r*T)*profit;
        discrete_data.push_back(discounted_profit);
    }

    return discrete_data;
}
//flag = true --> call
//flag = false --> put
void get_data_blocking_results(bool flag, Random& rnd, double ndata, double data_per_block, double nblocks, double sigma, double T, double nsteps, double S0, double K, double r){
    //Calcolo delle media di blocco
    vector<double>direct_data = get_direct(flag, rnd, ndata, sigma, T, S0, K, r);
    vector<double>discrete_data = get_discrete(flag, rnd, ndata, sigma, T, nsteps, S0, K, r);

    vector<double>block_means_dir, block_means_disc;
    block_means_dir.resize(nblocks);
    block_means_disc.resize(nblocks);

    for(int i = 0; i < nblocks; i++) {
        double sum_dir = 0.0;
        double sum_disc = 0.0;
        for(int j = 0; j < data_per_block; ++j) {
            sum_dir += direct_data[i*data_per_block + j];
            sum_disc += discrete_data[i*data_per_block + j];
        }
        block_means_dir[i] = sum_dir / data_per_block;
        block_means_disc[i] = sum_disc / data_per_block;
    }
    //Calcolo delle media progressive
    double cum_mean_dir  = 0.0;
    double cum_mean2_dir = 0.0;
    double cum_mean_disc  = 0.0;
    double cum_mean2_disc = 0.0;

    ofstream results_direct,results_discrete;

    string filename_direct = (flag)
                        ? "OUTPUT/results_direct_call.dat"
                        : "OUTPUT/results_direct_put.dat";
    string filename_discrete = (flag)
                        ? "OUTPUT/results_discrete_call.dat"
                        : "OUTPUT/results_discrete_put.dat";

    results_direct.open(filename_direct);
    results_discrete.open(filename_discrete);

    for(int b = 0; b < nblocks; ++b) {
        cum_mean_dir  += block_means_dir[b];
        cum_mean2_dir += block_means_dir[b] * block_means_dir[b];
        cum_mean_disc  += block_means_disc[b];
        cum_mean2_disc += block_means_disc[b] * block_means_disc[b];

        // numero di blocchi finora
        int m = b + 1;

        // media progressiva fino al blocco b
        double mean_b_dir  = cum_mean_dir  / m;
        double mean2_b_dir = cum_mean2_dir / m;
        double mean_b_disc  = cum_mean_disc  / m;
        double mean2_b_disc = cum_mean2_disc / m;

        double err_dir = (b == 0)
                        ? 0.0
                        : sqrt((mean2_b_dir - mean_b_dir * mean_b_dir) / (m - 1));
        double err_disc = (b == 0)
                        ? 0.0
                        : sqrt((mean2_b_disc - mean_b_disc * mean_b_disc) / (m - 1));

        results_direct << m << " " << mean_b_dir << " " << err_dir <<endl;
        results_discrete << m << " " << mean_b_disc << " " << err_disc <<endl;
        
    }
    results_discrete.close();
    results_direct.close();
}

int main() {
    Random rnd;
    int p1, p2;
    ReadPrimes("Primes", p1, p2);                     
    InitializeRandom(rnd, "seed.out", p1, p2);

    const int nsteps = 100; //Passi per discretizzazione
    const double S0 = 100.; //Valore iniziale del sottostante
    const double K = 100.; //Strike
    const double T = 1.; //Tempo di consegna
    const double r = 0.1; //Tasso di interesse senza rischio
    const double sigma = 0.25; //volatilità

    const int ndata = 1000000; //Simulazioni totali
    const int nblocks = 100; //Numero di blocchi
    const int data_per_block = ndata/nblocks; //Dati per blocco

    get_data_blocking_results(true, rnd, ndata, data_per_block, nblocks, sigma, T, nsteps, S0, K, r);
    get_data_blocking_results(false, rnd, ndata, data_per_block, nblocks, sigma, T, nsteps, S0, K, r);

   rnd.SaveSeed();
   return 0;
}