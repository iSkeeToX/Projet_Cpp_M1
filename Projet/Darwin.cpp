#include "Darwin.h"
#include "Matrix.h"

#include <iostream>
#include <random>
#include <cmath>
#include <fstream>
#include <filesystem> //Check if a file exist in the current directory
std::default_random_engine  re(time(0));

//_____________________________Fonctions

double recompense(const Matrix& MeanParameters){
    double Size = MeanParameters(0,0);
    double SizeHoles = MeanParameters(0,1);
    double Vol = MeanParameters(0, 2);
    double porosity = MeanParameters(0, 3);
    //double Surf_To_Vol_Ratio = MeanParameters(0, 4);
    double Sphericity = MeanParameters(0, 5);
    
    return -10*pow(Size - 6, 6) - 30*pow(SizeHoles - 1, 6) - pow(Vol - 7, 4) - 30*pow(porosity - 1./7, 4) - 10*pow(std::abs(Sphericity - sqrt(M_PI*sqrt(3))), 3);
}

//Fitness_Sortie_Du_Cul -10*pow(Size - 6, 6) - 10*pow(SizeHoles - 1, 6) - pow(Vol - 7, 4) - 10*pow(porosity - 1./7, 2) - 10*pow(std::abs(Sphericity - sqrt(M_PI*sqrt(3))), 3);
//Sigmoide_Du_Cul 1 / (1 + exp(10*pow(Size - 6, 6) + 10*pow(SizeHoles - 1, 6) + pow(Vol - 7, 4) + 10*pow(porosity - 1./7, 2) + 10*pow(std::abs(Sphericity - sqrt(M_PI*sqrt(3))), 3)));
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

int Pivot(Matrix& Scores, int low, int high, bool Ascending, int column){
    double pivot = Scores(high, column);
    int i = low - 1;
    if (Ascending){
        for(int k = low; k < high; k++){
            if(Scores(k, column) <= pivot){
                i++;
                std::swap(Scores(i, 0), Scores(k, 0));
                std::swap(Scores(i, 1), Scores(k, 1));
            }
        }
    }
    else{
        for(int k = low; k < high; k++){
            if(Scores(k, column) >= pivot){
                i++;
                std::swap(Scores(i, 0), Scores(k, 0));
                std::swap(Scores(i, 1), Scores(k, 1));
            }
        }
    }
    std::swap(Scores(i+1, 0), Scores(high, 0));
    std::swap(Scores(i+1, 1), Scores(high, 1));

    return i+1;
}

//Ascending = true -> Trie dans l'odre croissant
//colunm -> Trie selon la colonne column
void quicksort(Matrix& Scores, int low, int high, bool Ascending, int column){
    if(low < high){
        int pivot = Pivot(Scores, low, high, Ascending, column);

        quicksort(Scores, low, pivot - 1, Ascending, column);
        quicksort(Scores, pivot + 1, high, Ascending, column);
    }
}

void Darwin::Sort_Scores(){
    quicksort(Scores, 0, pop-1, false, 0);
}

void Darwin::Data(const std::string Name, double mean_crea, double stdev_crea, double acceptation) const{
    if (!(std::filesystem::exists(Name))){
        std::ofstream fich(Name);
        fich << "mean;stdev;acceptation;E(fitness);sigma(fitness);Generation;fitness1;Genes(Scores(0,1),:);fitness50;Genes(Scores(49,1),:);fitness100;Genes(Scores(99,1),:)\n";
        fich.close();
    }
    std::ofstream fich(Name, std::ios_base::app);
    
    fich << mean_crea << "; " << stdev_crea << "; " << acceptation << "; ";
    
    double mean = 0;
    for(int i=0; i<Scores.nx; i++){
        mean+=Scores(i,0);
    }
    mean = mean/Scores.nx;

    
    double stdev = 0;
    for(int i=0; i<Scores.nx; i++){
        stdev+= pow(Scores(i,0) - mean, 2);
    }
    stdev = sqrt(stdev/Scores.nx);

    fich << mean << "; " << stdev << "; " << generation << "; " << Scores(0,0) << "; [";
    
    int index = Scores(0,1);
    for(int i=0; i<Genes.ny; i++){
        fich << Genes(index, i) << ", ";
    }
    fich << "]; " << Scores(pop/2 - 1, 0) << "; ";

    index = Scores(pop/2 - 1, 1);
    for(int i=0; i<NbrGenes; i++){
        fich << Genes(index, i) << ", ";
    }
    fich << "]; " << Scores(pop - 1, 0) << "; ";

    index = Scores(pop - 1, 1);
    for(int i=0; i<NbrGenes; i++){
        fich << Genes(index, i) << ", ";
    }
    fich << "];\n";
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
            if(acc < acceptation){
                futurs_parents(i,0)=Scores(k,1);
            }
            else{ 
                futurs_parents(i,0)=Scores(j,1);
            }
        }
        else{ 
            if(acc < acceptation){
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
void Darwin::Nouveaux_genes(Matrix& futurs_parents, Matrix& couples, double mean, double stdev){
    int i,j,k;

    Matrix Genes_enfants = Matrix(2*couples.nx,21);

    Matrix enfant1 = Matrix(21,1);
    Matrix enfant2 = Matrix(21,1);

    std::normal_distribution<double> Normal(mean, stdev);
    std::uniform_int_distribution<int> Choix_Heredite(0, 1);
    std::uniform_real_distribution<double> Tirer_Proba_Mutation(0, 1);

    int shift = 0, choix_parent;
    for(i=0;i<couples.nx;i++){        
        for(j=0;j<21;j++){
            choix_parent = Choix_Heredite(re);

            if (Tirer_Proba_Mutation(re) < 0.001){
                enfant1(j, 0) = Normal(re);
            }
            else if (Tirer_Proba_Mutation(re) < 0.001){
                enfant1(j, 0) = -Genes(couples(i, choix_parent), j);
            }
            else{
                enfant1(j, 0) = Genes(couples(i, choix_parent), j);
            }

            if (Tirer_Proba_Mutation(re) < 0.001){
                enfant2(j, 0) = Normal(re);
            }
            else if (Tirer_Proba_Mutation(re) < 0.001){
                enfant2(j, 0) = -Genes(couples(i, 1-choix_parent), j);
            }
            else{
                enfant2(j, 0) = Genes(couples(i, 1-choix_parent), j);
            }
        }

        for(k=0;k<21;k++){
            Genes_enfants(shift,k)=enfant1(k,0);
            Genes_enfants(shift+1,k)=enfant2(k,0);
        }
        shift+=2;
    }

    quicksort(futurs_parents, 0, futurs_parents.nx - 1, true, 0);
    quicksort(Scores, 0, Scores.nx - 1, true, 1);
    k = 0;
    for(j=0; j < 21; j++){
        Genes(k, j) = Genes(futurs_parents(0 ,0), j);
    }
    Scores(k, 0) = Scores(futurs_parents(0 ,0), 0);
    Scores(k, 1) = k;
    k++;

    //Individus vainqueurs de Hunger Games
    for(i=1; i<futurs_parents.nx; i++){
        if (futurs_parents(i, 0) != futurs_parents(i-1, 0)){
            for(j=0; j < 21; j++){
                Genes(k, j) = Genes(futurs_parents(i, 0), j);
            }
            Scores(k, 0) = Scores(futurs_parents(i ,0), 0);
            Scores(k, 1) = k;
            k++;
        }
    }
    
    //Individus créés par les vainqueurs
    for(i = 0; i < Genes_enfants.nx; i++){
        for(j=0; j<21; j++){    
            Genes(k + i, j) = Genes_enfants(i, j);
        }
        Scores(k + i, 0) = 0;
    }

    //Potentiel reste d'individus si l'on a tiré 2 fois le même parent dans le tournoi
    for(i = k + Genes_enfants.nx; i<pop; i++){
        for(int j=0; j < 21; j++){
            Genes(i, j) = Normal(re);
        }
        Scores(i, 1) = i;
        Scores(i, 0) = 0;
    }
}

void Darwin::Next_Generation(double mean, double stdev, int Nparts, int N_Temps, int N_Steps, int N_Stat, std::string Name, double acceptation){
    
    for(int Individu=0; Individu < pop; Individu++){
        if(Scores(Individu, 0) == 0){
            (*this).Simulation(Individu, Nparts, N_Temps, N_Steps, N_Stat);
        }
        if(Individu == pop/2){
            std::cout << "L'individu " << Individu << " est simulé\n";
        }
        else if (Individu == 3*(pop/4)){
            std::cout << "L'individu " << Individu << " est simulé\n";
        }
    }
    std::cout << "Simulation de la génération " << generation << " terminée\n";
    (*this).Sort_Scores();

    (*this).Data(Name, mean, stdev, acceptation);
    
    Matrix Futurs_parents = (*this).futurs_parents(acceptation);
    Matrix Couples = (*this).couples(Futurs_parents);

    (*this).Nouveaux_genes(Futurs_parents, Couples, mean, stdev);
    std::cout << "La génération " << generation << " a enfanté la prochaine génération\n";
    generation++;
}


//_____________________________Le Constructeur de Darwin

Darwin::Darwin(int pop, int NbrGenes, double mean, double std, int nx, int ny) : pop(pop), generation(0), NbrGenes(NbrGenes), Genes(Matrix(pop, NbrGenes, mean, std)), Scores(Matrix(pop,2)), Ising(IsingModel(nx, ny)){

    for(int i=0; i < Scores.nx; i++){
        Scores(i,1) = i;
    }

}

Darwin::~Darwin(){}