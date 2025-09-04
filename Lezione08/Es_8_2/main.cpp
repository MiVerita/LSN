#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
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
//ritorna il valore della psi di test -- due gaussiane simmetriche
double get_psi(double x, double mu, double sigma){
   return exp(-pow(x-mu,2)/(2*sigma*sigma)) + exp(-pow(x+mu,2)/(2*sigma*sigma));
}
//ritorna il termine della derivata seconda della psi di test -- due gaussiane centrate
double get_second_derivative_psi(double x, double mu, double sigma){
   double t1 = exp(-pow(x-mu,2)/(2*sigma*sigma));
   double t2 = exp(-pow(x+mu,2)/(2*sigma*sigma));

   double d1 = (pow(x-mu,2)/pow(sigma,4) - 1/pow(sigma,2)) * t1;
   double d2 = (pow(x+mu,2)/pow(sigma,4) - 1/pow(sigma,2)) * t2;

   return (d1 + d2)/(t1 + t2);
}
//ritorna il valore di energia associata alla psi test -- due gaussiane centrate
double get_energy(double x, double mu, double sigma){
   return -0.5*get_second_derivative_psi(x, mu, sigma) + pow(x,4) - 2.5*x*x;
}
//container of VMC results
struct VMC_data{
   double energy; //Miglior valore medio di energia nello step VMC
   double error; //Errore associato al milgior valore medio di energia nello step VMC
   double accept_rate; //Acceptance in fase di misura nello step VMC
   vector<double> samples; //Valori campionati all'interno del passo VMC
   double final_delta; //Larghezza dello step finale all'interno del passo VMC
   double sigma; //Valore istantaneo di sigma
   double mu; //Valore istantaneo di mu
   vector<double> data_blocking_mean; //Medie progressive data blocking
   vector<double> data_blocking_error; //Erorri sulla media progressiva data blocking
};
//Esegue equilibrazione del parametro delta in modo da avere un'acceptance del 50%
double tune_delta(Random& rnd, double mu, double sigma, double input_delta){
   int equilibration = 5000; //Numero di passi di equilibrazione
   int check_acceptance = 100; //Equilibration/check_acceptance è il numero di controllo di delta
   double acceptance_interval = 0; //Acceptance per fase di controllo
   double target_acceptance = 0.5; //Acceptance voluta

   double x = 0.0; //Partenza algoritmo di metropolis
   double delta = input_delta;

   for (int i = 0; i < equilibration; i++) {
      double x_new = x + delta*(2.0*rnd.Rannyu() - 1.0); //proposta uniforme per il passo del metropolis
      double A = min(1.0, pow(get_psi(x_new,mu,sigma)/get_psi(x,mu,sigma),2));
      if (rnd.Rannyu() < A) {
         x = x_new;        // passo accettato
         acceptance_interval++;
      }

      if((i+1) % check_acceptance == 0){
         if((acceptance_interval/check_acceptance) < target_acceptance){
            delta *= 0.9; //se l'acceptance di blocco è più bassa di quella target diminuisco del 10% delta
         }else{
            delta*=1.1;  //se l'acceptance di blocco è più alta di quella target aumento del 10% delta
         }
         delta = max(0.1,min(delta, 5.0)); //limito i possibili valori di delta ta 0.1 e 5.0
         acceptance_interval = 0;
      }
   }

   return delta;
}
//Crea un campione di dati distribuiti come la densità data dalla psi di test
vector<double> get_samples(Random& rnd, double mu, double sigma, double delta, int nsteps){
   double acceptance_meas = 0.; //Acceptance in fase di misura
   double x = 0.0;
   vector<double> samples;
   //Generazione dei valori di energia tramite algoritmo di metropolis
   for(int i = 0; i < nsteps; i++){
      double x_new = x + delta*(2.0*rnd.Rannyu() - 1.0); //proposta uniforme per il passo del metropolis
      double A = min(1.0, pow(get_psi(x_new,mu,sigma)/get_psi(x,mu,sigma),2));
      if (rnd.Rannyu() < A) {
         x = x_new;        // passo accettato
         acceptance_meas++;
      }
      samples.push_back(x);
   }
   samples.push_back(acceptance_meas/nsteps);
   return samples;
}
/*
   Input:
      1. rnd -> oggetto random
      2. sigma -> valore istantaneo di sigma
      3. mu -> valore istantaneo di mu
      4. nsteps -> numero di passi di simulazione
      5. input_delta -> valore di delta iniziale, prima della fase di equilibrazione dell'algoritmo di metropolis
*/
//Esegue uno step di simulated annealing 
VMC_data VMC_step(Random& rnd, double sigma, double mu, int nsteps, double input_delta){
   double delta = tune_delta(rnd, mu, sigma, input_delta);

   //FASE DI MISURA
   vector<double> samples = get_samples(rnd, mu, sigma, delta, nsteps);
   double acceptance_meas = samples[nsteps];
   samples.pop_back();

   vector<double> energy;
   for(int i = 0; i < nsteps; i++){
      energy.push_back(get_energy(samples[i], mu, sigma));
   }

   //Statistica sui valori di energia
   int nblocks        = 1000;
   int data_per_block = nsteps / nblocks;

   // Calcolo delle medie di blocco istantanee
   vector<double> block_means(nblocks);
   for(int i = 0; i < nblocks; i++) {
      double sum = 0.0;
      for(int j = 0; j < data_per_block; ++j) {
         sum += energy[i*data_per_block + j];
      }
      block_means[i] = sum / data_per_block;
   }

   // Calcolo le media progressive delle medie di blocco ed i loro relativi errori
   vector<double> prog_mean(nblocks), prog_err(nblocks);
   double cum_mean  = 0.0;
   double cum_mean2 = 0.0;

   for(int b = 0; b < nblocks; ++b) {
      // accumulo somma e somma dei quadrati
      cum_mean  += block_means[b];
      cum_mean2 += block_means[b] * block_means[b];

      // numero di blocchi finora
      int m = b + 1;

      // media progressiva fino al blocco b
      double mean_b  = cum_mean  / m;
      double mean2_b = cum_mean2 / m;

      prog_mean[b] = mean_b;
      // errore su mean_b: varianza delle block_means divisa per (m-1)
      prog_err[b]  = (b == 0)
                     ? 0.0
                     : sqrt((mean2_b - mean_b * mean_b) / (m - 1));
   }
   VMC_data results = {prog_mean[nblocks-1], prog_err[nblocks-1], acceptance_meas, samples, delta, sigma, mu, prog_mean, prog_err};

   return results;
}

/*
   Input:
      1. rnd -> oggetto Random
      2. nsteps -> numero di "raffredamenti"
      3. first_sigma, first_mu -> parametri funzione test di partenza
      4. T0 -> temperatura iniziale
      5. cool -> fattore di raffredamento
      6. vmc_steps -> numero di passi di simulazione per raffredamento
      7. delta -> delta di input
*/

void simulated_annealing(Random& rnd, int ncycles, double first_sigma, double first_mu, double T0, double cool, int nsteps, double delta){
   ofstream sa_file("../OUTPUT/sa_data.dat");
   double sigma = first_sigma;            // Larghezze della funzione test
   double mu = first_mu;                  // centri della funzione test
   double T = T0;                         // Temperatura iniziale
   double best_E = 0;
   vector<double> best_means, best_errors; // Data blocking associato alla migliore stima di energia
   VMC_data best_data;

   for(int step=0; step < ncycles; step++) {
      // Calcola energia corrente nel punto (σ,μ)
      VMC_data current = VMC_step(rnd, sigma, mu, nsteps, delta);
      if(best_E > current.energy){
         best_data = current;
         best_E = current.energy;
      }
      // Scelgo nuovi valori di σ e μ perturbando quello vecchio con un valore gaussiano con μ_gauss = 0, σ_gauss = 0.1
      double sigma_new = sigma + rnd.Gauss(0.,0.1);
      double mu_new    = mu    + rnd.Gauss(0.,0.1);
      // Mantiene σ,μ > 0.1 per vincolare l'algoritmo
         if(sigma_new <= 0.1 || mu_new <= 0.1) { 
            sigma_new = sigma;
            mu_new    = mu;
         }

      VMC_data proposal = VMC_step(rnd, sigma_new, mu_new, nsteps, delta);
      double dE = proposal.energy - current.energy;
      // Accetto se migliora o con probabilità exp(-dE/T)
         if(dE < 0 || exp(-dE/T) > rnd.Rannyu()) {
            sigma = sigma_new;
            mu    = mu_new;
         }

      // Registra dati: passo, energia, errore, σ, μ, accettanza
      sa_file << step << " "
         << current.energy << " "
         << current.error  << " "
         << sigma          << " "
         << mu             << " "
         << current.accept_rate << endl;

      // Raffredda la temperatura
      T *= cool;
   }
   sa_file << "###" <<endl;
   sa_file << best_E << " " << best_data.error << " " << best_data.mu << " " << best_data.sigma <<endl;
   sa_file.close();
   cout << best_E << " " << best_data.error << " " << best_data.mu << " " << best_data.sigma <<endl;
   ofstream best_blocking, best_samples;

   best_blocking.open("../OUTPUT/best_blocking.dat");
   for(int i = 0; i < 100; i++){
      best_blocking << i << " " << best_data.data_blocking_mean[i] << " " << best_data.data_blocking_error[i] <<endl;
   }
   best_blocking.close();

   best_samples.open("../OUTPUT/best_parameters_samples.dat");
   for(unsigned int i = 0; i < best_data.samples.size(); i++){
      best_samples << best_data.samples[i] <<endl;
   }
   best_samples.close();

}


int main(){
   Random rnd;
   int p1, p2;
   ReadPrimes("Primes", p1, p2);                     
   InitializeRandom(rnd, "seed.out", p1, p2);
   
   double sigma = 1.; //sigma iniziale
   double mu    = 2.; //Mu iniziale  
   double delta = 2.; //Passo delta iniziale
   double ncycles = 200; //numero di passi SA
   double T0  = 6.; //Temperatura iniziale
   int nsteps = 100000; //Grandezza campione per ogni passo SA

   if(nsteps % 100 != 0){
      cerr << "PROBLEM: nsteps must be multiple of 100" << endl;
      return 0;
   }
   simulated_annealing(rnd, ncycles, sigma, mu, T0, 0.95, nsteps,delta);
   
   rnd.SaveSeed();
   return 0;
}