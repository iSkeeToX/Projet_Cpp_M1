#include "Darwin.h"
#include "Matrix.h"

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

//_____________________________Fonctions

double recompense(Matrix MeanParameters){
    MeanParameters(1,1);
    return 1.5;
}

//_____________________________Méthodes

void Darwin::Simulation(int Individu, int Nparts, int N_Temp, int N_Steps, int N_Stats){

    Ising.Initialise_Lattice(Nparts);
    double T_0 = std::abs(Genes(Individu, 0));

    int k=0;
    for(int j=0; j < 6; j++){
        for(int i=j; i< 6; i++){
            Ising.InteractionMap(i, j) = Genes(Individu, k);
            
            if( std::abs(Genes(Individu, k)) > T_0){
                T_0 = std::abs(Genes(Individu, k));
            }
            k++;
        }
    }

    Ising.beta = 1/T_0;
    Ising.Annealing(N_Temp, N_Steps);

    for(int stat=0; stat < N_Stats; stat++){
        Ising.Metropolis_Step();
    }

    Scores(Individu, 0) = recompense(ConComp(Ising).ClustersParameters().mean_columns());
    Scores(Individu, 1) = Individu;

    Ising.Vider_Lattice();
}

int Pivot(Matrix& Scores, int low, int high){
    double pivot = Scores(high, 0);
    int i = low - 1;

    for(int k = low; k < high; k++){
        if(Scores(k, 0) >= pivot){
            i++;
            std::swap(Scores(i, 0), Scores(k, 0));
            std::swap(Scores(i, 1), Scores(k, 1));
        }
    }
    std::swap(Scores(i+1, 0), Scores(high, 0));
    std::swap(Scores(i+1, 1), Scores(high, 1));

    return i+1;
}

void quicksort(Matrix& Scores, int low, int high){
    if(low < high){
        int pivot = Pivot(Scores, low, high);

        quicksort(Scores, low, pivot - 1);
        quicksort(Scores, pivot + 1, high);
    }
}

void Darwin::Sort_Scores(){
    quicksort(Scores, 0, pop-1);
}

void Darwin::Data(const std::string Name) const{
    std::ofstream fich(Name);
}

//_____________________________Le Constructeur de Darwin

Darwin::Darwin(int pop, int NbrGenes, double mean, double std, int nx, int ny) : pop(pop), generation(0), NbrGenes(NbrGenes), Genes(Matrix(pop, NbrGenes, mean, std)), Scores(Matrix(pop,2)), Ising(IsingModel(nx, ny)){

    for(int i=0; i < Scores.nx; i++){
        Scores(i,1) = i;
    }

}

Darwin::~Darwin(){}