#include <cmath>
#include <cstdlib>
#include <string>
#include <stdexcept>
#include <fstream>
#include <algorithm>
#include "System.h"

using namespace std;

System :: System(int first_n_genes, int first_n_chromosomes){
   n_genes = first_n_genes;
   n_chromosomes = first_n_chromosomes;
}
System :: ~System(){}

void System :: initialize_rnd(int world_rank){
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
           for(int i=0; i<4; ++i) {
             seed[i] += world_rank*53;  //modifico il seed in funzione del rank del processo in modo da avere generatori diversi
           }
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.out" << endl;
}
double System :: get_random(){
   return rnd.Rannyu();
}
double System :: get_random(double min, double max){
   return rnd.Rannyu(min, max);
}
void System::initialize_cities() {
   vector<city> storage; 
   city put;
   ifstream read;
   read.open("cap_prov_ita.dat");
   if (!read.is_open()) {
      cerr << "Errore: impossibile aprire cap_prov_ita.dat" << endl;
      exit(1);
   }
   double read_x, read_y;
   int iter = 1;
   while(read >> read_x >> read_y){
      put.x = read_x;
      put.y = read_y;
      put.id = iter;
      storage.push_back(put);
      iter++;
   }
   read.close();

   if ((int)storage.size() != n_genes) {
       cerr << "Attenzione: lette " << storage.size() << " città, ma n_genes = " << n_genes << endl;
   }

   set_cities(storage);
}
void System :: initialize_chromosomes(){
   vector<city> solution; //particular chromosome
   vector<city> city_dictionary; //vector with possible city to extract

   for(int k = 0; k < get_n_chromosomes(); k++){
      solution.clear();
      city_dictionary = get_cities();

      for(int i = 0; i < get_n_genes(); i++) {
         int idx = static_cast<int>(rnd.Rannyu(0., city_dictionary.size()));
         city choice = city_dictionary[idx];

         solution.push_back(choice);
         city_dictionary.erase(city_dictionary.begin() + idx);
      }
      chromosomes chromo;
      chromo.genes = solution;
      chromo.fitness = compute_L2(solution);

      population.push_back(chromo);
   }
}
void System :: finalize(){
   rnd.SaveSeed();
}

city System :: get_city_by_id(int id){
   for(city iter : cities){
      if(iter.id == id){
         return iter;
      }
   }
   cerr << "ERROR: city with ID " << id << " not found!" << endl;
   throw runtime_error("Invalid city ID");
}
void System :: set_n_genes(int new_n_genes){
   n_genes = new_n_genes;
}
void System :: set_n_chromosomes(int new_n_chromosome){
   n_chromosomes = new_n_chromosome;
}
void System :: set_cities(vector<city> new_cities){
   cities = new_cities;
}
void System :: set_population(vector<chromosomes> new_pop){
   population = new_pop;
}
int System :: get_n_genes(){
   return n_genes;
}
int System :: get_n_chromosomes(){
   return n_chromosomes;
}
vector<chromosomes> System :: get_population(){
   return population;
}
vector<city> System :: get_cities(){
   return cities;
}
double System :: get_best_half_mean(){
   vector<chromosomes> ordered = get_population();
   double mean = 0;
   if((n_chromosomes % 2) == 1 ){
      for(unsigned int i = 0; i < (n_chromosomes-1)/2; i++){
         mean += ordered[i].fitness;
      }
      mean /= (n_chromosomes-1)/2;
   }else{
      for(unsigned int i = 0; i < n_chromosomes/2; i++){
         mean += ordered[i].fitness;
      }
      mean /= n_chromosomes/2;
   }   
   return mean;
}
bool System :: check(chromosomes chromo){
   vector<city> genes = chromo.genes;
   vector<int> ids;
   for(int i = 0; i < get_n_genes(); i++){
      ids.push_back(genes[i].id);
   }
   sort(ids.begin(), ids.end());
   for(int i = 0 ; i < get_n_genes(); i++){
      if(ids[i] != i+1 ){
         return false;
      }
   }
   return true;
}
double System :: compute_L2(vector<city> solution){
   double L2 = 0.;
   vector<city> path = solution;
   for(int i = 0; i < get_n_genes(); i++){
      if(i == get_n_genes() - 1){
         L2 += pow(path[i].x - path[0].x,2) + pow(path[i].y - path[0].y,2);
      }else{
         L2 += pow(path[i].x - path[i+1].x,2) + pow(path[i].y - path[i+1].y,2);
      }
   }
   return L2;
}
vector<chromosomes> System :: order_pop(vector<chromosomes> unordered){
   vector<chromosomes> ordered = unordered;
   for (unsigned int i = 0; i < ordered.size(); i++) {
       int min_idx = i;
       for (unsigned int j = i + 1; j < ordered.size(); j++) {
           if (ordered[j].fitness < ordered[min_idx].fitness) {
               min_idx = j;
           }
       }
       if (min_idx != i) {
           swap(ordered[i], ordered[min_idx]);
       }
   }
   return ordered;
}

vector<city> System :: pair_permutation_mut(vector<city> solution){
   vector<city> mutation = solution;

   int idx = static_cast<int>(rnd.Rannyu(1., mutation.size()));
   if(idx == mutation.size() - 1){
      swap(mutation[idx], mutation[1]);  
   }else{
      swap(mutation[idx], mutation[idx+1]);
   }
   return mutation;
}
vector<city> System:: m_permutation_mut(vector<city> solution){
   int m = 1 + int(rnd.Rannyu() * ((n_genes - 1)/2 - 1)); // block dimension
   int start1 = 1 + int(rnd.Rannyu() * (n_genes - 2*m - 1)); // first block
   int start2 = start1 + m + int(rnd.Rannyu() * (n_genes - start1 - 2*m)); // second block
   
   for (int i = 0; i < m; ++i) {
       swap(solution[start1 + i], solution[start2 + i]);
   }
   return solution;
}
vector<city> System:: n_shift_mut(vector<city> solution) {
   int m = 1 + int(rnd.Rannyu() * (n_genes - 2)); // block dimension
   int n_start = 1 + int(rnd.Rannyu() * (n_genes - m - 1)); // starting point of the shift
   int n = 1 + int(rnd.Rannyu() * (n_genes - n_start - m - 1)); // shift 
   
   vector<city> block(solution.begin() + n_start, solution.begin() + n_start + m);
   solution.erase(solution.begin() + n_start, solution.begin() + n_start + m);
   solution.insert(solution.begin() + n_start + n, block.begin(), block.end());

   return solution;
}
vector<city> System :: inversion_mut(vector<city> solution){
   vector<city> mutation;
   mutation.push_back(solution[0]);
   for(unsigned int i = 0; i < solution.size() - 1; i++){
      mutation.push_back(solution[solution.size() - i - 1]);
   }
   return mutation;
}

chromosomes System :: selection(vector<chromosomes> current_pop){
      /*int p = 1.35;
      int size = static_cast<int>(current_pop.size());
   
      //int idx = static_cast<int>(size * pow(rnd.Rannyu(), p)) + 1;
      int idx = int(size * pow(rnd.Rannyu(), p));
      if (idx >= size) idx = size - 1; 
      chromosomes choice = current_pop[idx];
   
      return choice;
      */
    const int k = 3;  // dimensione del torneo
    int best_idx = static_cast<int>(rnd.Rannyu(0., current_pop.size()));
    
    // sfida k−1 “avversari” casuali
    for(int i = 1; i < k; ++i) {
        int challenger = static_cast<int>(rnd.Rannyu(0., current_pop.size()));
        if (current_pop[challenger].fitness < current_pop[best_idx].fitness) {
            best_idx = challenger;
        }
    }
    
    return current_pop[best_idx];

}
vector<chromosomes> System :: crossover(chromosomes ex_parent_1, chromosomes ex_parent_2){
   chromosomes offspring_1, offspring_2;
   vector<chromosomes> offsprings;

   vector<city> save_1;
   vector<city> save_2;

   chromosomes parent_1 = ex_parent_1;
   chromosomes parent_2 = ex_parent_2;

   int n_division = static_cast<int>(rnd.Rannyu(1.,n_genes - 1));

   //saves transfer material
   for(unsigned int i = 0; i < n_genes - n_division; i++){
      save_1.push_back(parent_1.genes[n_division+i]);
      save_2.push_back(parent_2.genes[n_division+i]);
   }
   for(unsigned int i = 0; i < n_genes; i++){
      if(i >= n_division){
         auto contains = [](const vector<city>& genes, const city& c) { //lambda expression to control if a gene already exist in the chromosome
            for(const auto& g : genes){
               if(g.id == c.id) return true;
            }
            return false;
         };

         // Completamento di offspring_1 con i geni mancanti dal partner (parent_2)
         for(unsigned int j = 0; j < parent_2.genes.size(); j++){
            if(!contains(offspring_1.genes, parent_2.genes[j])){
               offspring_1.genes.push_back(parent_2.genes[j]);
            }
         }

         // Completamento di offspring_2 con i geni mancanti da parent_1
         for(unsigned int j = 0; j < parent_1.genes.size(); j++){
            if(!contains(offspring_2.genes, parent_1.genes[j])){
               offspring_2.genes.push_back(parent_1.genes[j]);
            }
         }

      }else{
         offspring_1.genes.push_back(parent_1.genes[i]);
         offspring_2.genes.push_back(parent_2.genes[i]);
      }
   }
   offsprings.push_back(offspring_1);
   offsprings.push_back(offspring_2);
   return offsprings;
}

void System::generation() {
   vector<chromosomes> old_generation = order_pop(get_population());
   vector<chromosomes> ordered_population = old_generation;
   vector<chromosomes> new_generation, new_gen_cross;

   // Elitism - saves top 2%
   int elite_count = static_cast<int>(n_chromosomes / 50);
   new_generation.insert(new_generation.end(), old_generation.begin(), old_generation.begin() + elite_count);
   old_generation.erase(old_generation.begin(), old_generation.begin() + elite_count);

   // Crossover
   int size = (old_generation.size() % 2 == 0) ? old_generation.size() : old_generation.size() - 1;

   for (int i = 0; i < size / 2; ++i) {
       chromosomes parent_1 = selection(ordered_population);
       chromosomes parent_2 = selection(ordered_population);

       //50% chance which selected parents get crossed
       if (rnd.Rannyu() < 0.7) {
           vector<chromosomes> offsprings = crossover(parent_1, parent_2);
           new_gen_cross.push_back(offsprings[0]);
           new_gen_cross.push_back(offsprings[1]);
       } else {
           new_gen_cross.push_back(parent_1);
           new_gen_cross.push_back(parent_2);
       }
   }

   // if dimension were odd put the last chromosome
   if (old_generation.size() % 2 == 1) {
       new_gen_cross.push_back(old_generation.back());
   }

   // Mutation -- each non elite can get mutated with 10% chance. Multiple mutation may occur
   for (unsigned int i = 0; i < new_gen_cross.size(); i++) {
       if (rnd.Rannyu() < 0.5)
           new_gen_cross[i].genes = pair_permutation_mut(new_gen_cross[i].genes);
       if (rnd.Rannyu() < 0.5)
           new_gen_cross[i].genes = n_shift_mut(new_gen_cross[i].genes);
       if (rnd.Rannyu() < 0.5)
           new_gen_cross[i].genes = m_permutation_mut(new_gen_cross[i].genes);
       if (rnd.Rannyu() < 0.5)
           new_gen_cross[i].genes = inversion_mut(new_gen_cross[i].genes);

       new_gen_cross[i].fitness = compute_L2(new_gen_cross[i].genes);
   }

   // Combine elite with new chromosomes
   new_generation.insert(new_generation.end(), new_gen_cross.begin(), new_gen_cross.end());
   new_generation = order_pop(new_generation);
   // refresh population
   population = new_generation;
}

vector<double> System :: serialize(vector<chromosomes> chromo){
   vector<double> serialize;
   for(unsigned int i = 0; i < chromo.size(); i++){
      int size = chromo[i].genes.size();

      for(unsigned int j = 0; j < size; j++){
         serialize.push_back(chromo[i].genes[j].id);
      }
      serialize.push_back(chromo[i].fitness);
   }   
   return serialize;
}
vector<chromosomes> System :: build_back(vector<double> serialized){
   vector<chromosomes> chromo_back;
   int n_chromosomes = serialized.size() / (n_genes + 1);
   for(unsigned int i = 0; i < n_chromosomes; i++){
      chromosomes single_chromo;
      for(int j = 0; j < n_genes; j++){
         single_chromo.genes.push_back(get_city_by_id(serialized[i*(n_genes+1) + j]));
      }
      single_chromo.fitness = serialized[i * (n_genes + 1) + n_genes];
      chromo_back.push_back(single_chromo);
   }
   return chromo_back;
}
void System :: refresh_best_pop(vector<chromosomes> best_pop){
   vector<chromosomes> pop = get_population(); 
   pop.erase(pop.end()-best_pop.size(), pop.end());
   pop.insert(pop.end(), best_pop.begin(), best_pop.end());
   set_population(pop);
}
