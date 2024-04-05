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

void EvolutionDeDarwin(int NbrGenes, double mean, double stdev){
    int pop =10, nx = 30, ny = 30, Nparts = (nx*ny)/9;
    int N_Temps = 100, N_Steps = nx*ny, N_stats = nx*ny;//Aucune de ces valeurs n'est utilisée normalement, c'est principalement pour quele programme aille vite
    int acceptation = 0.9;
    
    Darwin D = Darwin(MutationGaussiannFlip, PlusDeTrous, pop, NbrGenes, mean, stdev, nx, ny);

    for(int Individu = 0; Individu < pop; Individu++){
        D.Simulation(Individu, Nparts, N_Temps, N_Steps, N_stats);
    }
    cout << "Le comportement de tous les individus est simulé\n";
    
    D.Sort_Scores();
    
    //D.data(Name, mean, stdev, acceptation) Ecrit les données dans un fichier texte

    cout << "Gènes de la première génération :\n" << D.Genes;

    Matrix futurs_parents = D.futurs_parents(acceptation);
    Matrix Couples = D.couples(futurs_parents);

    D.Nouveaux_genes(futurs_parents, Couples, mean, stdev);

    quicksort(futurs_parents, 0 , futurs_parents.nx -1, true, 0);
    cout << "Voici les futur_parents :\n" << futurs_parents;
    cout << "Les gènes de la nouvelle génération :\n" << D.Genes;

    cout << "Et voilà, nous sommes passés d'une génération à l'autre\n";
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
        
        cout << "Fitness : " << PremierTest(ConComp(Ising).ClustersParameters().mean_columns()) << "\n";
        cout << "Fitness Simul : " << fitness << '\n';
        Ising.TakePicture();
    }
    fich.close();
}

//Ceci est le main que nous utilisions pour lancer des simulations
void maint(){

    int nx = 30, ny = 30;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    int pop = 100, NbrGenes = 21;

    double mean = 0, stdev = 10;

    double acceptation = 0.9;


    cout << "mean : " << mean << "stdev : " << stdev << "\n";
    Darwin D = Darwin(MutationGaussiannFlip, TestGazExp, pop, NbrGenes, mean, stdev, nx, ny);

    for(int i=0; i<100; i++){
        D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "TestGazExpGaussien.txt", acceptation);
    }

    if(false){
    Darwin D = Darwin(MutationFlip, TestGaz, pop, NbrGenes, nx, ny);
    double mean = 0, stdev = 0;
    for(int i=0; i<100; i++){
        D.Next_Generation(mean, stdev, Nparts, N_Temps, N_Steps, N_Stat, "TestGaz.txt", acceptation);
    }
    }
}

int main(){
    int nx = 50, ny = 50;
    int Nparts = (nx*ny)/9;

    int N_Temps = 100;
    int N_Steps = 10 * Nparts;

    int N_Stat = 100 * Nparts;

    int NbrGenes = 21;
    double mean = 0, stdev = 10; //kT

    float taille = 10;

    cout << "Bonjour, voici le projet de Le Bon Louison et Jouve Alexandre intitulé :\n";
    cout << "Principes d'auto assemblage pour des particules avec des géométries simple et des interactions complexes\n";
    cout << "Ce travail se base sur la thèse éponyme de Lara Koehler.\n";

    cout << "Ceci est un menu (si si c'en est un) afin de vous donner une vision un peu éclatée du projet.\n";
    cout << "Les petites fonctions que vous allez faire tourner ici sont écrites dans le fichier main.cpp\n\n";

    //cout << "0 : Notre note au projet C++\n",
    cout << "1 : Une illustration de notre algorithme de Métropolis\n";
    cout << "2 : Une démonstration des méthodes principales de la classe ConComp\n";
    cout << "3 : Une démonstration de la vie d'une génération\n\n";

    int ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage = 0;
    
    cout << "Merci de renseigner un nombre entier compris entre 1 et 3\n";
    cin >> ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage;

    if(ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage == 1){
        cout << "Tirons un individu au hasard selon une loi normale de paramètre N(0 kT,10 kT)\n";
        cout << "Vous allez assiter à la simulation complète de cet individu\n";
        cout << "Cette simulation comprend une phase d'Annealing et une phase de statistique\n";
        Metropolis(nx, ny, Nparts, taille, 10*N_Temps, N_Steps, N_Stat);
    }
    else if (ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage == 2){
        cout << "Nous allons remplir aléatoirement un Lattice de taille 50*50\n";
        cout << "Isoler ses composantes connexes et ensuite Identifier puis isoler la plus grande composante connexe\n";
        cout << "Vous pouvez fermer chaque fenêtre qui apparaît devant vous lorsque vous le souhaitez, cela passera immédiatement à l'étape suivante\n";   
        ConnectedComponentDisplay(nx, ny, Nparts, taille);
    }
    else if(ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage == 3){
        cout << "Il va se dérouler sous vos yeux la magnifique entreprise de la vie, l'évolution de Darwin\n";
        cout << "Nous allons simuler une génération d'une population tirée aléatoirement, essayez de voir les nouveaux nés ainsi que les parents conservés\n";
        EvolutionDeDarwin(NbrGenes, mean, stdev);
    }

    cout << "Souhaitez vous regarder autre chose ? (0 ou 1)\n";
    cout << "0 : Non\n";
    cout << "1 : Oui\n";

    cin >> ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage;

    if(ChoixAffuteEtRationnelDeMKebailiMaisUnPeuEmotionnelIlDoitDecrypter150ProjetEn4JoursQuelCourage != 1){
        return 0;
    }
    else{
        cout << "Au revoir, et à bientôt pour la présentation orale !\n";
        main();
    }
    //AfficherPretendants("GazzExp.txt", nx, ny, Nparts, N_Temps, N_Steps, N_Stat);
    //ConnectedComponentDisplay(nx, ny, Nparts, taille);
    //50*50, 277 particules, Carte d'interaction Gaussienne N(0, 10)
    //Metropolis(nx, ny, Nparts, taille, N_Temps, N_Steps, N_Stat);


    //Characterize_Fluct(false, false);

}

