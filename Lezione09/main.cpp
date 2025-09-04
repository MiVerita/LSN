#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "System.h"


using namespace std;

int main(){
   ofstream results_cities;
   ofstream results_genetic;
   ofstream results_best_ch;

   const int n_cities = 34;
   const int n_chromosomes = 200;
   bool flag = true; //true -- circle <> false -- square

   if(flag){
      results_cities.open("OUTPUT/results_cities_circle");

      if(!results_cities.is_open()){
         cerr << "PROBLEM: Unable to open results_cities_circle" << endl;
         return 0;
      }
      results_genetic.open("OUTPUT/results_genetic_circle");
      
      if(!results_genetic.is_open()){
         cerr << "PROBLEM: Unable to open results_genetic_circle" << endl;
         return 0;
      }
      results_best_ch.open("OUTPUT/results_best_ch_circle");

      if(!results_best_ch.is_open()){
         cerr << "PROBLEM: Unable to open results_best_ch_circle" << endl;
         return 0;
      }

   }else{
      results_cities.open("OUTPUT/results_cities_square");

      if(!results_cities.is_open()){
         cerr << "PROBLEM: Unable to open results_cities_square" << endl;
         return 0;
      }
      results_genetic.open("OUTPUT/results_genetic_square");
      
      if(!results_genetic.is_open()){
         cerr << "PROBLEM: Unable to open results_genetic_square" << endl;
         return 0;
      }

      results_best_ch.open("OUTPUT/results_best_ch_square");

      if(!results_best_ch.is_open()){
         cerr << "PROBLEM: Unable to open results_best_ch_square" << endl;
         return 0;
      }
   }

   System vitro(n_cities, n_chromosomes);
   
   vitro.initialize_rnd();
   vitro.initialize_cities(flag);
   vitro.initialize_chromosomes();
   
   for(unsigned int i = 1; i < 9000; i++){
      vitro.generation();
      results_genetic << i << " " << vitro.get_best_half_mean() <<endl;
      results_best_ch << i << " " << vitro.get_population()[0].fitness <<endl;
   }

   vector<chromosomes> population = vitro.get_population();
   for(city iter : population[0].genes){
      results_cities << iter.x << " " << iter.y <<endl;
   }
   results_cities << population[0].genes[0].x << " " << population[0].genes[0].y <<endl;
   
   vitro.finalize();
   results_genetic.close();
   results_cities.close();
   results_best_ch.close();

   return 0;
}