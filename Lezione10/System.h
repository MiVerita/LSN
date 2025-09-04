#ifndef __System__
#define __System__

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include "random.h"

using namespace std;

struct city{
    int id; //id della città
    double x;   //posizione x della città
    double y;   //posizione y della città
};
struct chromosomes{
    double fitness; //valore della loss per la soluzione
    vector<city> genes; //array di geni rappresentanti una proposta di soluzione
};

class System{
    private: 
        int n_genes; //chromosomes size
        int n_chromosomes; //population size
        vector<chromosomes> population; //population storage
        vector<city> cities; //cities dictionary
        Random rnd;

    public:
        System(int first_n_genes, int first_n_chromosomes);
        ~System();

        //initialization methods
        void initialize_rnd(int world_rank);
        void initialize_cities();
        void initialize_chromosomes(); //initialize the first population -- flag = false --> square; flag = true --> circumference
        void finalize();

        //set/get methods
        void set_n_genes(int new_n_genes); //set n_genes value
        void set_n_chromosomes(int new_n_chromosome); //set n_chromosome value
        void set_cities(vector<city> storage); //set cities vector
        void set_population(vector<chromosomes> new_pop); //set population vector
        int get_n_genes();  //return n_genes value
        int get_n_chromosomes(); //return n_chromosome value
        vector<chromosomes> get_population(); //get population vector
        vector<city> get_cities(); //get cities vector
        double get_best_half_mean();
        city get_city_by_id(int id);
        double get_random();
        double get_random(double min, double max);
        
        //operations on chromosomes methods
        bool check(chromosomes chromo); //check if a chromosome is valid
        vector<chromosomes> order_pop(vector<chromosomes> unordered); //modify population's vector by sorting per fitness quality: between two chromosomes the one with lower value of L2 is better
        double compute_L2(vector<city> solution); //computes the L1 fitness
        chromosomes selection(vector<chromosomes> pop); //select a chromosome based on the fitness attribute of the chromosomes
        vector<chromosomes> crossover(chromosomes parent_1, chromosomes parent_2); //gives two offspring from two chromosomes selected
        vector<city> pair_permutation_mut(vector<city> solution); //mutates a chromosome by pair permutation of ids in the solution
        vector<city> m_permutation_mut(vector<city> solution); //mutates a chromosome by m (casual number) permutation of ids in the solution
        vector<city> n_shift_mut(vector<city> solution); //mutates a chromosome by n shift of ids in the solution
        vector<city> inversion_mut(vector<city> solution); //mutates a chromosome by inverting the ids in the solution
        
        vector<double> serialize(vector<chromosomes> chromo); //serialize a chromosomes to be sent via communicator
        vector<chromosomes> build_back(vector<double> chromo_back); //build back a serialized chromosome
        void refresh_best_pop(vector<chromosomes> best_pop); //update the curent continent population with wich passed via communicator
        
        void generation(); //computes a generation
};
#endif