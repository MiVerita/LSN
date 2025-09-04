#include <iostream>
#include <fstream>
#include <mpi.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "System.h"

using namespace std;

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  int world_size, world_rank;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  // 1) Ogni processo costruisce il proprio Continente
  bool flag = false;
  const int n_chromosomes = 250;
  const int n_cities = 110;
  const int Ngen = 10000;
  const int Nmigr = 500;       // ogni 500 generazioni migrano

  System vitro(n_cities, n_chromosomes);
  vitro.initialize_cities();
  vitro.initialize_rnd(world_rank); //Ogni rank inizializza il proprio generatore in funzione del rank chiamante in modo da avere generatori diversi per ogni rank
  vitro.initialize_chromosomes();
  
  // 2) Ciclo evolutivo con migrazioni
  ofstream cout_best_L2, cout_best_L2_mean, cout_best_path;
  if(flag){
    cout_best_L2.open("OUTPUT/best_L2_migration.dat", ios::app);
    cout_best_L2_mean.open("OUTPUT/best_L2_mean_migration.dat", ios::app);
    cout_best_path.open("OUTPUT/best_path_migration.dat", ios::app);
  }else{
    cout_best_L2.open("OUTPUT/best_L2_indipendent.dat", ios::app);
    cout_best_L2_mean.open("OUTPUT/best_L2_mean_indipendent.dat", ios::app);
    cout_best_path.open("OUTPUT/best_path_indipendent.dat", ios::app);
  }
  //INIZIO GENERAZIONI
  for(int gen = 1; gen <= Ngen; gen++){

    if(gen % 100 == 0){
      cout << "RANK: " << world_rank << "; GENERAZIONE: " << gen <<endl;
    }

    vitro.generation();
    if(flag){
      if(gen % Nmigr == 0){
        //GENERAZIONE DEI PARTNERS DI SCAMBIO
        vector<int> partners(world_size);
        //IL RANK 0 E' IL RANK MASTER -- genera i partner
        if (world_rank == 0) {
            vector<int> ranks(world_size);
            for (int i = 0; i < world_size; ++i)
                ranks[i] = i;
        
            // Genera un vettore mescolato dei rank
            vector<int> shuffle_ranks;
            for (int i = 0; i < world_size; ++i) {
                int choice = static_cast<int>(vitro.get_random(0., ranks.size()));
                shuffle_ranks.push_back(ranks[choice]);
                ranks.erase(ranks.begin() + choice);
            }
        
            // Crea accoppiamenti simmetrici
            for (int i = 0; i < world_size; i += 2) {
                partners[shuffle_ranks[i]]     = shuffle_ranks[i + 1];
                partners[shuffle_ranks[i + 1]] = shuffle_ranks[i];
            }
        }

        //Il rank master 0 informa gli altri rank dei partner di scambio scelti
        MPI_Bcast(partners.data(), partners.size(), MPI_INT, 0, MPI_COMM_WORLD); 

        //La migrazione scambia il 10% migliore di ogni popolazione
        vector<chromosomes> migrants;
        for(unsigned int i = 0; i < 5; i++){
          migrants.push_back(vitro.get_population()[i]);
        }

        vector<double> serialized_migrants = vitro.serialize(migrants);
        vector<double> serialized_answer(serialized_migrants.size());

          //INVIO DEI MESSAGGI DAI VARI RANK; RICEZIONE MESSAGGI DAI VARI RANK
        MPI_Sendrecv(
          serialized_migrants.data(), serialized_migrants.size(), MPI_DOUBLE, partners[world_rank], 0,
          serialized_answer.data(), serialized_migrants.size(), MPI_DOUBLE, partners[world_rank], 0,
          MPI_COMM_WORLD, MPI_STATUS_IGNORE
        );
          //RICOSTRUZIONE MESSAGGIO RICEVUTI DAI RANK
        vector<chromosomes> answer;
        answer = vitro.build_back(serialized_answer);
        vitro.refresh_best_pop(answer);
        vitro.set_population(vitro.order_pop(vitro.get_population()));

        cout << "RANK: " << world_rank << "; MIGRAZIONE DELLA GENERAZIONE: " << gen << " ESEGUITA DA " << world_rank << " A " << partners[world_rank] <<endl;
      }
    }
    double local_best_L2 = vitro.get_population()[0].fitness;
    double local_L2_mean = vitro.get_best_half_mean();
    double global_best_L2;
    double global_best_L2_mean;
    
    MPI_Reduce(&local_best_L2, &global_best_L2, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_L2_mean, &global_best_L2_mean, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);

    //IL RANK 0 E' IL RANK MASTER -- stampa i risultati
    if(world_rank == 0){
      cout_best_L2 << gen << " " << global_best_L2 << endl;
      cout_best_L2_mean << gen << " " << global_best_L2_mean << endl;
    }

    if(gen == Ngen){
      vector<chromosomes> local_best_chromo;
      local_best_chromo.push_back(vitro.get_population()[0]);

      vector<double> local_best_chromo_serialized = vitro.serialize(local_best_chromo);
      vector<double> gathered_chromo((n_cities+1)*world_size);
      vector<chromosomes> global_best_chromo;

      //IL RANK 0 E' IL RANK MASTER -- ottiene il miglior percorso
      MPI_Gather(local_best_chromo_serialized.data(), local_best_chromo_serialized.size(), MPI_DOUBLE,
      gathered_chromo.data(), local_best_chromo_serialized.size(), MPI_DOUBLE,
        0, MPI_COMM_WORLD);

      if(world_rank == 0){
        vector<double> get_fitness;

        for(int i = n_cities; i < gathered_chromo.size(); i += n_cities+1){
          get_fitness.push_back(gathered_chromo[i]);
        }

        double min_fit = get_fitness[0];
        int index_min = 0;

        // Trovo il miglior fitness e a quale cromosoma appartiene
        for (int i = 1; i < get_fitness.size(); ++i) {
            if (get_fitness[i] < min_fit) {
              min_fit = get_fitness[i];
              index_min = i;
            }
        }
        //Ricostruisco il cromosoma migliore
        for (int i = world_size - 1; i >= 0; --i) {
            if (i != index_min) {
                int start = i * (n_cities + 1);
                int end = start + (n_cities + 1);
                gathered_chromo.erase(gathered_chromo.begin() + start, gathered_chromo.begin() + end);
            }
        }
        global_best_chromo = vitro.build_back(gathered_chromo);

        for(chromosomes iter_chromo : global_best_chromo){
          for(city iter_gene : iter_chromo.genes){
            cout_best_path << iter_gene.x << " " << iter_gene.y <<endl;
          }
        }

      }
    }
  }

  vitro.finalize();

  cout_best_L2.close();
  cout_best_L2_mean.close();
  cout_best_path.close();

  MPI_Finalize();
  return 0;
}