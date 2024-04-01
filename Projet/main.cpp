#include "Matrix.h"
#include "TriangularLattice.h"
#include "IsingModel.h"


#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <algorithm> //min,max
#include <vector>
#include <fstream>

using namespace std;

#include<random>

void FaitTournerMetropolis() {
    
    int nx=256, ny=256;
    int Nparts=(nx*ny)/9;
    float Taille=8;

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

    Ising.Annealing(N_T, N_steps, Taille, window);
    
    Ising.Metropolis_Step();
}

void TesterConComp(){
    int nx=20, ny=20;
    float Taille=10;

    IsingModel L=IsingModel(nx, ny);
    L.Initialise_Lattice(100);

    ConComp Connected = ConComp(L);

    Connected.write("Connected_Components.txt");
    sf::Color custom(127, 127, 127, 255); //Gris 

    Matrix Size = Connected.SizeConComps();
    
    for(int ConCompNumber=1; ConCompNumber < Connected.NbrCC; ConCompNumber++){
        cout << "ConComp : " << ConCompNumber << " Length : " << Connected.OuterBorderLength(ConCompNumber) << " Size : " << Size(ConCompNumber-1,0) << "\n";
    }

    Connected.Show_Connected_Components(2*Taille);
    cout << Connected.NbrCC << "\n";
}

int main(){
    ConComp Connected = ConComp(100, 100);
    Connected.NbrCC = 1;

    int x = 10, y = 10;
    Connected[Connected.site_xy(x,y)]=1, Connected[Connected.site_xy(x+1,y)]=1, Connected[Connected.site_xy(x+2,y)]=1;
    Connected[Connected.site_xy(x+1,y+1)]=1, Connected[Connected.site_xy(x+1,y+2)]=1;
    Connected[Connected.site_xy(x+2,y+2)]=1, Connected[Connected.site_xy(x+3,y+2)]=1,Connected[Connected.site_xy(x+3,y+1)]=1, Connected[Connected.site_xy(x+3,y)]=1;
    
    int xp = x + 4;
    Connected[Connected.site_xy(xp,y)]=1, Connected[Connected.site_xy(xp+1,y)]=1, Connected[Connected.site_xy(xp+2,y)]=1;
    Connected[Connected.site_xy(xp+1,y+1)]=1, Connected[Connected.site_xy(xp+1,y+2)]=1;
    Connected[Connected.site_xy(xp+2,y+2)]=1, Connected[Connected.site_xy(xp+3,y+2)]=1,Connected[Connected.site_xy(xp+3,y+1)]=1, Connected[Connected.site_xy(xp+3,y)]=1;
    
    Matrix Parameters = Connected.ClustersParameters();

    cout << Parameters;

    Connected.Show_Connected_Components(2*10);
    
}
