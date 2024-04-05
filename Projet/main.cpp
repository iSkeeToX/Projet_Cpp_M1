#include "Matrix.h"
#include "TriangularLattice.h"
#include "IsingModel.h"
#include "Darwin.h"

#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <algorithm> //min,max
#include <vector>
#include <fstream>

using namespace std;

#include<random>

//Cette fonction illustre une simulation typique d'un individu avec l'algorithme de Monte-Carlo
void Metropolis(int nx, int ny, int Nparts, float taille, int N_T, int N_Steps, int N_Stats){

    IsingModel Ising = IsingModel(nx, ny);

    Ising.Gaussian_InteractionMap(0, 10);

    cout << "Carte d'interaction : " << Ising.InteractionMap;

    Ising.Initialise_Lattice(Nparts);

    sf::Color custom(127, 127, 127, 255); //Gris 
    sf::RenderWindow Window(sf::VideoMode(880, 760), "Illustration Algorithme de Monte-Carlo (Annealing)");
    
    Ising.Annealing(N_T, N_Steps, taille, Window);

    int step = 0;

    sf::RenderWindow window(sf::VideoMode(880,760), "Illustration Algorithme de Monte-Carlo");
    while((window.isOpen())){
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        window.clear(custom);
        for(int i = 0; (i < pow(10, 4)) && (step < N_Stats); i++){
            Ising.Metropolis_Step();
            step++;
        }
        Ising.affiche_SFML(window, taille);
        window.display();
    }
}

//Cette fonction illustre les fonctionnalités principales de la classe ConComp
//Vous verrez agir le constructeur de la classe, qui passe d'un Lattice à ses composantes connexes
//Puis nous isolerons la composante connexe de plus grande taille et nous calculerons sont périmètre
void ConnectedComponentDisplay(int nx, int ny, int Nparts, float taille){

    IsingModel Ising= IsingModel(nx, ny);
    Ising.Initialise_Lattice(Nparts);

    sf::Color custom(127, 127, 127, 255); //Gris 
    sf::RenderWindow window(sf::VideoMode(880,760), "Illustration Classe ConComp");
    
    while((window.isOpen())){
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear(custom);
        Ising.affiche_SFML(window, taille);
        window.display();
    }

    ConComp Connected = ConComp(Ising);
    

    Matrix Size = Connected.SizeConComps();
    int CCN = 0, SizeN = 0;
    for(int i=0; i<Connected.NbrCC; i++){
        if(Size(i, 0) > SizeN){
            CCN = i + 1;
            SizeN = Size(i, 0);
        }
    }

    ConComp Isolated = Connected.isolateConComp(CCN);
    std::cout << "Nous allons désormais isoler la plus grande composante connexe qui se trouve être la composante connexe numéro " << CCN << "\n";
    std::cout << "La longueur du périmètre de cette composante connexe est : " << Isolated.OuterBorderLength(CCN) << "\n";
    Connected.Show_Connected_Components(1.2*taille);
    Isolated.Show_Connected_Components(3*taille);
}

//Cette fonction nous permet d'effectuer les captures d'écran des Individus que l'on a sélectionné manuellement lors d'une simulation
void AfficherPretendants(string Name, int nx, int ny, int Nparts, int N_Temps, int N_Steps, int N_Stat){
    ifstream fich(Name);

    if (!fich.is_open()) {
        std::cerr << "Error opening file: " << Name << std::endl;
        return;
    }
    double mean, stdev, acceptation, fitness, meanfitness;
    
    fich >> mean >> stdev >> acceptation;
     
    for(int k = 0; k<3; k++){
        IsingModel Ising = IsingModel(nx, ny);
        fich >> meanfitness >> fitness;
        
        Ising.Initialise_Lattice(Nparts);

        for(int j=0; j<6; j++){
            for(int i=j; i<6; i++){
                fich >> Ising.InteractionMap(i, j);
            }
        }

        Ising.Annealing(N_Temps, N_Steps);

        for(int k=0; k< N_Stat; k++){
            Ising.Metropolis_Step();
        }
        
        cout << "Fitness : " << recompense(ConComp(Ising).ClustersParameters().mean_columns()) << "\n";
        cout << "Fitness Simul : " << fitness << '\n';
        Ising.TakePicture();
    }
    
    fich.close();

}

//Cette fonction nous permet de Caractériser, pour une carte d'interaction donnée l'espérance et l'écart type
//de la fitness après simulation, cela nous permet d'effectuer les graphes situés dans /Data/Fluct/
void Characterize_Fluct(bool GaussianOrNot, bool SameConfigOrNot){
    
    std::ofstream fich("FluctFitness.dat");
    
    std::random_device rdd;
    std::mt19937 genn(rdd());

    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;
    
    double mean=0, stdev = 10;
    double T_0i = 10, T_0;

    IsingModel Ising=IsingModel(nx, ny);
    
    if(GaussianOrNot){//InteractionMap Gaussienne
    std::normal_distribution Gaussian(mean ,stdev);
    Ising.Gaussian_InteractionMap(mean, stdev);
    std::ofstream fich1("IntMap.dat");

    for(int j=0; j<6; j++){
        for(int i=j; i<6; i++){
            fich1 << Ising.InteractionMap(i, j) << " ";
            if( Ising.InteractionMap(i, j) > T_0i){
                T_0i = std::abs(Ising.InteractionMap(i, j));
            }
        }
    }
    fich1.close();
    }
    else{//InteractionMap Discrète dans -10, 0, 10;
    std::ofstream fich1("IntMapDiscrete.dat");
    std::uniform_int_distribution<int> unif(-1,1);

        for(int j=0; j<6; j++){
        for(int i=j; i<6; i++){
            Ising.InteractionMap(i, j) = 10*unif(rdd);
            fich1 << Ising.InteractionMap(i, j) << " ";
            if( Ising.InteractionMap(i, j) > T_0i){
                T_0i = std::abs(Ising.InteractionMap(i, j));
            }
        }
    }
    fich1.close();
    }

    if(SameConfigOrNot){//La configuration initiale change à chaque itération
    for(int i=0; i < 100; i++){
        T_0 = T_0i;
        Ising.Initialise_Lattice(Nparts);
        Ising.beta = 1/T_0;
        Ising.Annealing(N_Temps, N_Steps);

        for(int stat=0; stat < N_Stat; stat++){
            Ising.Metropolis_Step();
        }
        fich << recompense(ConComp(Ising).ClustersParameters().mean_columns()) << "\n";
        Ising.Vider_Lattice();
        fich.close();
    }
    }
    else{//La configuration initiale est la même pour tout le monde
    IsingModel Config = IsingModel(nx, ny);
    Config.Initialise_Lattice(Nparts);

    for(int i=0; i < 100; i++){
        T_0 = T_0i;
        Ising.Initialise_Lattice(Nparts);
        Ising = Config;
        Ising.Particles = Config.Particles;

        Ising.beta = 1/T_0;
        Ising.Annealing(N_Temps, N_Steps);

        for(int stat=0; stat < N_Stat; stat++){
            Ising.Metropolis_Step();
        }
        fich << recompense(ConComp(Ising).ClustersParameters().mean_columns()) << "\n";
        Ising.Vider_Lattice();
        fich.close();
    }
    }
}


void aa(){

    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    int pop = 100, NbrGenes = 21;

    //double mean = 0, stdev = 20;

    double acceptation = 0.9;

    if(false){
    double stdev = 10;
    for(double mean : {-10, 0, 10}){
        cout << "mean : " << mean << "stdev : " << stdev << "\n";
        Darwin D = Darwin(MutationGaussiannFlip, recompense, pop, NbrGenes, mean, stdev, nx, ny);

        for(int i=0; i<100; i++){
            D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "New_Aim.txt", acceptation);
        }
    }

    stdev = 8;
    for(double mean : {-5, 0, 5}){
        cout << "mean : " << mean << "stdev : " << stdev << "\n";
        Darwin D = Darwin(MutationGaussiannFlip, recompense, pop, NbrGenes, mean, stdev, nx, ny);

        for(int i=0; i<100; i++){
            D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "New_Aim.txt", acceptation);
        }
    }
    }

    if(false){
    Darwin D = Darwin(MutationFlip, recompense, pop, NbrGenes, nx, ny);
    double mean = 0, stdev = 0;
    for(int i=0; i<500; i++){
        D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "DiscreteFlipTrous4.txt", acceptation);
    }
    }
    //AfficherPretendants("test.txt", nx, ny, Nparts, N_Temps, N_Steps, N_Stat);

    if(true){
    Darwin D = Darwin(MutationFlip, TestGaz, pop, NbrGenes, nx, ny);
    double mean = 0, stdev = 0;
    for(int i=0; i<100; i++){
        D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "TestGaz.txt", acceptation);
    }
    }
}

int main(){
    if(true){
    
    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    float taille = 10;

    //AfficherPretendants("test.txt", nx, ny, Nparts, N_Temps, N_Steps, N_Stat);
    //ConnectedComponentDisplay(nx, ny, Nparts, taille);
    //50*50, 277 particules, Carte d'interaction Gaussienne N(0, 10)
    //Metropolis(nx, ny, Nparts, taille, N_Temps, N_Steps, N_Stat);
    }

    Characterize_Fluct(false, false);

}

