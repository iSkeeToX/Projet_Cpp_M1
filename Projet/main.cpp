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
void Metropolis(){
    int nx=50, ny=50;
    int Nparts = (nx*ny)/9;

    float taille = 10;

    IsingModel Ising = IsingModel(nx, ny);

    int N_T = 100; //Number of temperatures
    int N_Steps = 10 * (nx * ny); //Number of steps by temperature
    int N_Stats = 100 *(nx*ny); //Number of steps on which we "average"

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


void FaitTournerMetropolis() {
    
    int nx=256, ny=256;
    int Nparts=(nx*ny)/9;
    //float Taille=8;

    IsingModel Ising=IsingModel(nx, ny);
    Ising.Initialise_Lattice(Nparts);
    Ising.Gaussian_InteractionMap(0,10);

    Ising.InteractionMap.afficher();
    cout << "___________\n";

    int N_T = 100; //Number of temperatures
    int N_steps = pow(10,4); //Number of steps by temperature
    //int N_statistics=pow(10,4); //Number of steps on which we do the averaging


    sf::Color custom(127, 127, 127, 255); //Gris 
    sf::RenderWindow window(sf::VideoMode(1500,900), "Mon super projet");

    //Ising.beta = ...
    Ising.Annealing(N_T, N_steps);
    //Ising.Annealing(N_T, N_steps, Taille, window);
    
    Ising.Metropolis_Step();
}

void ConCompTest(){
    int nx=20, ny=20;
    float Taille=10;

    IsingModel L= IsingModel(nx, ny);
    L.Initialise_Lattice(100);

    ConComp Connected = ConComp(L);

    sf::Color custom(127, 127, 127, 255); //Gris 

    Matrix Size = Connected.SizeConComps();
    
    for(int ConCompNumber=1; ConCompNumber < Connected.NbrCC; ConCompNumber++){
        cout << "ConComp : " << ConCompNumber << " Length : " << Connected.OuterBorderLength(ConCompNumber) << " Size : " << Size(ConCompNumber-1,0) << "\n";
    }

    Connected.Show_Connected_Components(2*Taille);
    cout << Connected.NbrCC << "\n";
}

void COnCompTest(){
    IsingModel Ising = IsingModel(25,50);

    Ising.Initialise_Lattice(250);
    
    ConComp Connected = ConComp(Ising);
    //Connected.NbrCC = 1;

    //int x = 10, y = 10;
    //Connected[Connected.site_xy(x,y)]=1, Connected[Connected.site_xy(x+1,y)]=1, Connected[Connected.site_xy(x+2,y)]=1;
    //Connected[Connected.site_xy(x+1,y+1)]=1, Connected[Connected.site_xy(x+1,y+2)]=1;
    //Connected[Connected.site_xy(x+2,y+2)]=1, Connected[Connected.site_xy(x+3,y+2)]=1,Connected[Connected.site_xy(x+3,y+1)]=1, Connected[Connected.site_xy(x+3,y)]=1;
    
    //int xp = x + 4;
    //Connected[Connected.site_xy(xp,y)]=1, Connected[Connected.site_xy(xp+1,y)]=1, Connected[Connected.site_xy(xp+2,y)]=1;
    //Connected[Connected.site_xy(xp+1,y+1)]=1, Connected[Connected.site_xy(xp+1,y+2)]=1;
    //Connected[Connected.site_xy(xp+2,y+2)]=1, Connected[Connected.site_xy(xp+3,y+2)]=1,Connected[Connected.site_xy(xp+3,y+1)]=1, Connected[Connected.site_xy(xp+3,y)]=1;
    
    cout << Connected.NbrCC << "\n";

    Matrix Parameters = Connected.ClustersParameters();

    cout << Parameters;

    Connected.Show_Connected_Components(2*7);
    
}

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

void maint(){
    
    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    int pop = 100, NbrGenes = 21;

    //double mean = 0, stdev = 20;

    double acceptation = 0.9;

    double stdev = 10;
    for(double mean : {-10, 0, 10}){
        cout << "mean : " << mean << "stdev : " << stdev << "\n";
        Darwin D = Darwin(pop, NbrGenes, mean, stdev, nx, ny);

        for(int i=0; i<100; i++){
            D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "Plus_De_Trou.txt", acceptation);
        }
    }

    stdev = 8;
    for(double mean : {-5, 0, 5}){
        cout << "mean : " << mean << "stdev : " << stdev << "\n";
        Darwin D = Darwin(pop, NbrGenes, mean, stdev, nx, ny);

        for(int i=0; i<100; i++){
            D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "Plus_De_Trou.txt", acceptation);
        }
    }
    
    

    //AfficherPretendants("test.txt", nx, ny, Nparts, N_Temps, N_Steps, N_Stat);

}

void mmain(){
    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    AfficherPretendants("test.txt", nx, ny, Nparts, N_Temps, N_Steps, N_Stat);
}

int main(){
    //50*50, 277 particules, Carte d'interaction Gaussienne N(0, 10)
    //Metropolis();

    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    AfficherPretendants("test.txt", nx, ny, Nparts, N_Temps, N_Steps, N_Stat);
}