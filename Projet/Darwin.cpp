#include "Darwin.h"
#include "Matrix.h"

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>

std::default_random_engine  re(time(0));

//_____________________________Fonctions

double recompense(Matrix MeanParameters){
    MeanParameters(0,1);
    std::uniform_real_distribution<double> fitness(0,100);
    return fitness(re);
}

//_____________________________Méthodes

void Darwin::Simulation(int Individu, int Nparts, int N_Temp, int N_Steps, int N_Stats){

    Ising.Initialise_Lattice(Nparts);
    double T_0 = std::abs(Genes(Individu, 0));

    int k=0;
    for(int j=0; j < 6; j++){
        for(int i=j; i < 6; i++){
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

//Creation des couples de parents
Matrix Darwin::futurs_parents(const double acceptation) const{
    
    int i,j,k;
    double acc;

    
    Matrix futurs_parents = Matrix(pop/2, 2); //2 juste pour pas refaire l'algo de sort
    
    std::uniform_int_distribution<int> Indidivu_Aleatoire(0, pop - 1);
    std::uniform_real_distribution<double> Tirer_Proba_Acceptation(0, 1);
    
    for(i=0; i<pop/2; i++){   // tournoi
        k=Indidivu_Aleatoire(re);
        j=Indidivu_Aleatoire(re);    
        acc= Tirer_Proba_Acceptation(re);

        if (Scores(k,0)>=Scores(j,0)){
            if(acc > acceptation){
                futurs_parents(i,0)=Scores(k,1);
            }
            else{ 
                futurs_parents(i,0)=Scores(j,1);
            }
        }
        else{ 
            if(acc > acceptation){
                futurs_parents(i,0)=Scores(j,1);
            }
            else{ 
                futurs_parents(i,0)=Scores(k,1);
            }
        }
    }

    return futurs_parents;
}

Matrix Darwin::couples(const Matrix& futurs_parents) const{
    Matrix couples = Matrix(pop/4,2);
    std::uniform_int_distribution<int> Tirage_Parent(0,  pop/2 - 1);
    int k, j;

    for(int i=0; i< pop/4 ; i++){   // on considère négligeable la proba de tirer deux fois le même couple et même si c'était le cas pas forcéments mêmes enfants  
        
        k=Tirage_Parent(re);
        j=Tirage_Parent(re);

        couples(i,0)=futurs_parents(k,0);
        couples(i,1)=futurs_parents(j,0);
    
    }

    return couples;
}

//Creation de la nouvelle génération de gênes
void Darwin::Nouveaux_genes(Matrix& futurs_parents, Matrix& couples){
    int i,j,k;

    Matrix Genes_enfants = Matrix(2*couples.nx,21);
    Matrix mixing = Matrix(21,1);
    
    Matrix enfant1 = Matrix(21,1);
    Matrix enfant2 = Matrix(21,1);


    std::uniform_int_distribution<int> Choix_Heredite(0, 1);
    std::uniform_real_distribution<double> Tirer_Proba_Mutation(0, 1);

    for(i=0;i<couples.nx;i++){
        
        
        for(j=0;j<21;j++){
            mixing(j,0)=Choix_Heredite(re);
        }
       

        for(j=0;j<21;j++){
            if(mixing(j,0)==1){

                if (Tirer_Proba_Mutation(re) < 0.001 ){
                    enfant1(j,0)=-Genes(couples(i,0), j);
                }
                else{
                    enfant1(j,0)=Genes(couples(i,0), j);  
                }

                if (Tirer_Proba_Mutation(re) < 0.001){
                    enfant2(j,0)=-Genes(couples(i,1), j);
                }
                else{
                    enfant2(j,0)=Genes(couples(i,1), j);
                }
            }
            else{

                if (Tirer_Proba_Mutation(re) < 0.001){
                     enfant1(j,0)=-Genes(couples(i,1), j);
                }
                else{
                   enfant1(j,0)=Genes(couples(i,1), j);
                }
                if (Tirer_Proba_Mutation(re) < 0.001){
                    enfant2(j,0)=-Genes(couples(i,0),j);  
                }
                else{
                    enfant2(j,0)=Genes(couples(i,0),j);  
                }          
            }
        
            for(k=0;k<21;k++){
                Genes_enfants(i,k)=enfant1(k,0);
                Genes_enfants(i+1,k)=enfant2(k,0);
            }
        }
    }

    quicksort(futurs_parents, 0, futurs_parents.nx - 1);
    k = 0;
    for(j=0; j < 21; j++){
        Genes(k, j) = Genes(futurs_parents(futurs_parents.nx - 1,0), j);
    }
    k++;

    for(i=1; i<futurs_parents.nx; i++){
        if (futurs_parents(futurs_parents.nx - i, 0) != futurs_parents(futurs_parents.nx - 1 - i, 0)){
            for(j=0; j < 21; j++){
                Genes(k, j) = Genes(futurs_parents(futurs_parents.nx - i - 1,0), j);
            }
            k++;
        }
    }
    
    for(i = 0; i < Genes_enfants.nx; i++){
        for(j=0; j<21; j++){    
            Genes(k + i, j);
            //Genes_enfants(i, j);
        }
    }
}

void Darwin::Next_Generation(int Nparts, int N_Temps, int N_Steps, int N_Stat, std::string Name, double acceptation){
    
    for(int Individu=0; Individu < pop; Individu++){
        (*this).Simulation(Individu, Nparts, N_Temps, N_Steps, N_Stat);
        std::cout << "L'individu " << Individu << "est simulé\n";
    }
    (*this).Sort_Scores();

    (*this).Data(Name);
    

    Matrix Futurs_parents = (*this).futurs_parents(acceptation);
    Matrix Couples = (*this).couples(Futurs_parents);

    (*this).Nouveaux_genes(Futurs_parents, Couples);
    generation++;
}


//_____________________________Le Constructeur de Darwin

Darwin::Darwin(int pop, int NbrGenes, double mean, double std, int nx, int ny) : pop(pop), generation(0), NbrGenes(NbrGenes), Genes(Matrix(pop, NbrGenes, mean, std)), Scores(Matrix(pop,2)), Ising(IsingModel(nx, ny)){

    for(int i=0; i < Scores.nx; i++){
        Scores(i,1) = i;
    }

}

Darwin::~Darwin(){}