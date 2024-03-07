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

void Garder() {
    
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

int main(){
    int nx=100, ny=100;
    float Taille=10;

    IsingModel L=IsingModel(nx, ny);
    L.Initialise_Lattice(2500);
    ConComp Connected = ConComp(L);

    Connected.write("Connected_Components.txt");
    sf::Color custom(127, 127, 127, 255); //Gris 
    sf::RenderWindow window(sf::VideoMode(600,500), "Mon super projet");
    while(window.isOpen()){
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        window.clear(custom);
        L.affiche_SFML(window, Taille);
       
        window.display();
    }

    Connected.Show_Connected_Components(Taille);
    cout << Connected.NbrCC << "\n";
    
}
