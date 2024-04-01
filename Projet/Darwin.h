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
    
    Matrix couples(double proportion_parents, double acceptation ); //retourne la liste des positions des futurs parents 
    Matrix Nouveaux_genes(int nbr_enfants, double proportion_enfants,Matrix Couples);
    
    Darwin(int pop, int NbrGenes, double mean, double std, int nx, int ny);
    ~Darwin();


};











#endif