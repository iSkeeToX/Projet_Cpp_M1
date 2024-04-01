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









//Creation des couples de parents

Matrix Darwin::couples(double proportion_parents, double acceptation ){
    int i,j,k,l;
    Matrix(Darwin.pop*proportion_parents/2,2) couples;
    //Matrix(Darwin.pop/4,2) couples_deja_selectionnes;
    //for(i=0,i<Darwin.pop/4,i++){
        //couples_deja_selectionnes(i,0)=100;
        //couples_deja_selectionnes(i,1)=100;
    //}
    Matrix(Darwin.pop*proportion_parents, 1) futurs_parents;
    std::default_random_engine  re(time(0));
    for(i=0;i<Darwin.pop*proportion_parents;i++){   // tournoi
        k=std::uniform_int_distribution<int> distrib{  0,  Darwin.pop  };
        j=std::uniform_int_distribution<int> distrib{  0,  Darwin.pop  };    
         
        if (Scores(k,0)>=Scores(j,0))
            double acc= uniform_distribution<double> distrib{ 0, 1 };
            if(acc>acceptation)
                futurs_parents(i,0)=Scores(k,1);
            else 
                futurs_parents(i,0)=Scores(j,1);
        else 
            double acc= uniform_distribution<double> distrib{ 0, 1 };
            if(acc>acceptation)
                futurs_parents(i,0)=Scores(j,1);
            else 
                futurs_parents(i,0)=Scores(k,1);

            
        }
    for(i=0;i<Darwin.pop*proportion_parents/2;i++){   // on considère négligeable la proba de tirer deux fois le même couple et même si c'était le cas pas forcéments mêmes enfants  
        
        k=std::uniform_int_distribution<int> distrib{  0,  Darwin.pop*proportions_parents  };
        j=std::uniform_int_distribution<int> distrib{  0,  Darwin.pop*proportion_parents  };

        couples(i,0)=futurs_parents(k,0);
        couples(i,0)=futurs_parents(j,0);
    
    }

return couples;

    
    

}



//Creation de la nouvelle génération de gênes

void Darwin::Nouveaux_genes(int nbr_enfants, double proportion_enfants, Matrix couples){
    int i,j,k;
    Matrix(2*couples.nx,21) Genes_enfants;
    for(i=0;i<couples.nx;i++){
        Matrix(21,1) mixing;
        for(j=0;j<21;j++){
            mixing(j,0)=std::uniform_int_distribution<int> distrib{  0, 1  };
        }
        Matrix(21,1) enfant1;
        Matrix(21,1) enfant2;
        for(j=0;j<21;j++){
        if(mixing(j,0)==1){
            double mutation=uniform_distribution<double> distrib{  0, 1  }
            if (mutation<0.001){
                enfant1(j,0)=-Genes(couples(i,0),j);
            }
            else{
                enfant1(j,0)=Genes(couples(i,0),j);  
            }
            double mutation2=uniform_distribution<double> distrib{  0, 1  }
            if (mutation2<0.001){
                enfant2(j,0)=-Genes(couples(i,1),j);
            }
            else{
                enfant2(j,0)=Genes(couples(i,1),j);
            }
        
        else{
            double mutation3=uniform_distribution<double> distrib{  0, 1  }
            if (mutation3<0.001){
                enfant1(j,0)=-Genes(couples(i,1))
            }
            else{
                enfant1(j,0)=Genes(couples(i,1),j);
            }
            double mutation4=uniform_distribution<double> distrib{  0, 1  }
            if (mutation4<0.001){

                enfant2(j,0)=-Genes(couples(i,0),j);  
            }
            else{
                enfant2(j,0)=Genes(couples(i,0),j);  
            }          
        }
        }   
        for(k=0;k<21;k++){
            Genes_enfants(i,k)=enfant1(k,0);
            Genes_enfants(i+1,k)=enfant1(k,0);
        }
    }
    Matrix(pop*proportion_parents,0) indices_parents;
    for(i=0;i<pop*proportion_parents;i++){ 
        if (i%2==0){ 
        indices_parents(i,0)=couples(i,0);
        indices_parents(i+1,0)=couples(i,1);
        }
    for(i=0;i<pop;i++){
        int iter=0;
        for(j=0;j<pop*proportion_parents;j++){
            if(j==indices_parents(j,0)){
                continue
            }
            else{
                for(k=0;k<21;k++){
                    Genes(i,k)=Genes_enfants(iter,0);
                }
                iter+=1;
            }
        }



        }
    
        }   

    }

}
//_____________________________Le Constructeur de Darwin

Darwin::Darwin(int pop, int NbrGenes, double mean, double std, int nx, int ny) : pop(pop), generation(0), NbrGenes(NbrGenes), Genes(Matrix(pop, NbrGenes, mean, std)), Scores(Matrix(pop,2)), Ising(IsingModel(nx, ny)){

    for(int i=0; i < Scores.nx; i++){
        Scores(i,1) = i;
    }

}

Darwin::~Darwin(){}