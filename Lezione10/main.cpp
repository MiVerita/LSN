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
  const int Nmigr = 50;       // ogni 500 generazioni migrano

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

/*
  for(unsigned int j = 0; j < vitro.get_population().size(); j++){
    output << "RECEIVE -- > RANK: " << world_rank << "; ELEMENT: " << j << " ";
    for(unsigned int i = 0; i < vitro.get_population()[j].genes.size(); i++){
      output << vitro.get_population()[j].genes[i].id << " ";
    }
    output << "Fitness: " << vitro.get_population()[j].fitness <<endl;
  }
*/
12.516 38.0174
9.66781 39.9705
9.34238 40.1277
9.11331 39.2172
8.54034 39.3075
8.76505 39.4937
8.67964 40.0266
8.92197 40.7778
9.50674 40.9253
10.3403 42.7902
11.2878 42.7751
11.7639 43.5171
11.0941 43.9357
10.8687 43.9741
10.4544 44.0178
10.4019 43.7159
10.0526 44.2131
10.5941 44.6087
10.9924 45.4385
10.6708 45.1693
10.037 45.2209
10.3281 44.8014
9.66653 44.8476
9.49168 45.2613
9.75422 45.7567
9.18963 45.4642
9.27883 45.6395
9.13783 45.0369
8.74503 44.835
7.68249 45.0678
7.31966 45.7371
8.3836 46.1156
8.08691 45.567
7.55814 44.4581
8.25257 44.2334
7.86674 43.9584
8.93386 44.4073
8.20269 44.826
8.34628 45.5554
8.546 45.5843
8.75413 45.8397
9.14941 45.9395
9.41202 45.9005
9.69123 44.2384
10.936 44.5385
11.3426 44.4938
11.4064 45.6349
10.4259 45.7796
10.2584 46.3234
11.1258 46.0664
11.8279 44.7668
12.2742 44.9772
12.631 43.9465
12.0141 44.0227
12.059 44.3641
11.8734 45.4077
12.2063 45.8067
12.6597 45.9563
12.0789 46.2805
11.2302 46.6559
12.3346 45.4372
13.2358 46.0635
13.6252 45.9441
13.7773 45.6496
12.7014 43.6941
12.403 43.107
13.1509 43.153
13.2188 43.4801
13.6388 43.0922
13.5396 42.8834
13.9576 42.3103
13.6979 42.6581
12.8859 42.4147
12.4829 41.8933
12.4397 42.6539
11.4676 43.1672
11.2556 43.7698
11.9488 42.493
13.0126 41.4595
13.5758 41.6285
13.6103 42.1369
14.4159 42.1027
14.2081 41.6495
14.7057 41.2476
15.4529 41.5028
14.8262 41.7173
14.1169 41.2035
14.2488 40.8359
15.1406 40.9965
15.3106 40.4194
15.8216 40.5173
16.1466 41.1802
16.862 41.1258
17.0806 40.5488
17.6884 40.6359
18.2261 40.1522
16.8783 39.1874
16.4736 40.4476
16.3331 39.5967
16.4316 38.83
16.0987 38.6267
15.5542 38.1938
14.2807 37.5668
15.0874 37.5024
15.6398 38.1035
15.2907 37.0646
14.7213 36.922
14.0632 37.4899
13.5747 37.3123
13.3524 38.1112
12.516 38.0174