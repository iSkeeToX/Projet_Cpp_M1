#ifndef _DARWIN_H_
#define _DARWIN_H_

#include "Matrix.h"
#include "IsingModel.h"

class Darwin{
    public:
    
    int pop, generation, NbrGenes;
    Matrix Genes; //Genes(i,:) -> gènes du frérot i
    Matrix Scores; //Liste triée par fitness décroissante Genes(Scores(i,1),:) -> Gènes du ième meilleur

    //Fonction fitness, cette fonction est cruciale pour le bon fonctionnement de l'algorithme générationnel
    //Elle permet de noter les différents individus, nous avons choisi d'essayer de la maximiser
    //Malheureusement, il n'y a pas de méthode générale pour créer une fonction de récompense, nous devons alors tatonner
    double (*FitnessFunction)(const Matrix& MeanParameters);

    //Définit la manière dont mute les enfants lorsqu'ils naissent
    Matrix (*mutation)(const Matrix& Genes, const Matrix& couples, double mean, double stdev);



    IsingModel Ising;

    void Simulation(int Individu, int Nparts, int N_Temps, int N_Steps, int N_Stat);
    void Sort_Scores();
    void Data(const std::string Name, double mean, double stdev, double acceptation) const;
    
    Matrix futurs_parents(const double acceptation) const;
    Matrix couples(const Matrix& futurs_parents) const; //retourne la liste des positions des futurs parents 
    void Nouveaux_genes(Matrix& futurs_parents,  Matrix& couples, double mean, double stdev);
    
    void Next_Generation(double mean, double stdev, int Nparts, int N_Temps, int N_steps, int N_Stat, std::string Name, double acceptation);

    Darwin(Matrix mutation(const Matrix&, const Matrix&, double, double), double FitFunc(const Matrix&), int pop, int NbrGenes, double mean, double std, int nx, int ny);
    Darwin(Matrix mutation(const Matrix&, const Matrix&, double, double), double FitFunc(const Matrix&), int pop, int NbrGenes, int nx, int ny);

    ~Darwin();


};

double recompense(const Matrix& MeanParameters);

Matrix MutationGaussiannFlip(const Matrix& Genes, const Matrix& couples, double mean, double stdev);

Matrix MutationFlip(const Matrix& Genes, const Matrix& couples, double mean, double stdev);




#endif