#ifndef _DARWIN_H_
#define _DARWIN_H_

#include "Matrix.h"
#include "IsingModel.h"

class Darwin{
    public:
    
    int pop, generation, NbrGenes;
    Matrix Genes; //Genes(i,:) -> gènes du frérot i
    Matrix Scores; //Liste triée par fitness décroissante Genes(Scores(i,1),:) -> Gènes du ième meilleur

    IsingModel Ising;

    void Simulation(int Individu, int Nparts, int N_Temps, int N_Steps, int N_Stat);
    void Sort_Scores();
    void Data(const std::string Name) const;

    Darwin(int pop, int NbrGenes, double mean, double std, int nx, int ny);
    ~Darwin();


};











#endif