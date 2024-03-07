#include "TriangularLattice.h"
#include "Matrix.h"

#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <algorithm> //min,max
#include <vector>

using namespace std;


#include<random>



//Initialise le lattice L avec Nparticules ainsi qu'une matrice Nparts*2 contenant les coordonnées des particules (x,y)
Matrix Initialise_Lattice(Lattice& L, int Nparts){
    
    Matrix ParticlesCoordinates = Matrix(Nparts,2);
    
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()  

    std::uniform_int_distribution<int8_t> orientation(1, 6);
    std::uniform_int_distribution<int> SiteAlea(0, L.nx*L.ny-1);
    
    int SiteAl=SiteAlea(gen);
    int8_t Or=orientation(gen);

    for(int i=0;i<Nparts;i++){
        SiteAl=SiteAlea(gen);
        if(L[L.site_index(SiteAl)]==0){
            Or=orientation(gen);
            L[L.site_index(SiteAl)]=Or;

            ParticlesCoordinates(i,0)= SiteAl/L.ny;//Coord x
            ParticlesCoordinates(i,1)= (SiteAl + L.ny)%L.ny;
        }
        else{
            i--;
        }
    }
    
    return ParticlesCoordinates;
}

Matrix initialise_face_count(const Lattice& L){
    Matrix N_faces=Matrix(6,6);
    int site = 0;
    int face = 0;

    int min = 0;
    int max = 0;

    for(int i=0;i<L.nx;i++){
        for(int j=0;j<L.ny;j++){
            site = (int) L[L.site_xy(i,j)];
            if(site != 0){
                face = 0;
                for(Site voisin: L.voisins(L.site_xy(i,j))){
                    if((int) L[voisin] != 0){
                        max = std::max((site-1+6-face)%6,( ((int) L[voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
                        min = std::min((site-1+6-face)%6,( ((int) L[voisin])-1+3-face+6)%6);                       
                        N_faces(max,min) += 0.5;
                    }
                    face++;
                }
            }
        }
    }
    return N_faces;
}

Matrix Contact_Faces(const Lattice& L, const Site s, vector<int>& EmptySites){
    Matrix Faces = Matrix(6,2);
    int face=0;
    int max = 0,min = 0;

    int site = (int) L[s];

    for(Site voisin: L.voisins(s)){
        if( (int) L[voisin] == 0 ){
            Faces(face,0) = -1;
            Faces(face,1) = -1;
            EmptySites.push_back(face);
        }
        else{
            max = std::max((site-1+6-face)%6,( ((int) L[voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
            min = std::min((site-1+6-face)%6,( ((int) L[voisin])-1+3-face+6)%6);
            Faces(face,0)=max;
            Faces(face,1)=min;
        }
        face++;
    }
    return Faces;
}

Matrix Contact_Faces(const Lattice& L, const Site s){
    Matrix Faces = Matrix(6,2);
    int face=0;
    int max = 0,min = 0;

    int site = (int) L[s];

    for(Site voisin: L.voisins(s)){
        if( (int) L[voisin] == 0 ){
            Faces(face,0) = -1;
            Faces(face,1) = -1;
        }
        else{
            max = std::max((site-1+6-face)%6,( ((int) L[voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
            min = std::min((site-1+6-face)%6,( ((int) L[voisin])-1+3-face+6)%6);
            Faces(face,0)=max;
            Faces(face,1)=min;
        }
        face++;
    }
    return Faces;
}

void Metropolis_step(Lattice& L, Matrix& Particles, const Matrix InteractionMap, float beta){
        std::random_device rd;  // Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd() 


        std::uniform_int_distribution<int> RandomParticle(0, Particles.nx-1);
        int Particle_Number = RandomParticle(gen);
        Site s=L.site_xy(Particles(Particle_Number, 0), Particles(Particle_Number, 1));
    
        int Orientation = L[s];
        double DeltaE=0;

        //Calcul de la contribution du site s à l'énergie
        vector<int> EmptySites; //Liste contenant l'indice des sites vides (ds l'ordre trig)  
        Matrix FacesInitial = Contact_Faces(L,s, EmptySites);
        for(int i=0;i<6;i++){
            if (FacesInitial(i,0) >= 0){
                DeltaE-=InteractionMap(FacesInitial(i,0), FacesInitial(i,1));
            }
        }


        std::uniform_real_distribution<float> distribution(0, 1);
        float MoveOrSpin = distribution(gen);
    
        Site NewLocation = s;
        int NewOrientation = Orientation;


        if( (EmptySites.size() != 0) && (MoveOrSpin < 0.5) ){
            uniform_int_distribution<int> MoveWhere(0,EmptySites.size()-1);
            NewLocation = L.voisins(s)[EmptySites[MoveWhere(gen)]];
        
            L[s]=0;
            L[NewLocation] = Orientation;      
        }
        else{
            uniform_int_distribution<int> OrientationHow(1,5);
            NewOrientation = OrientationHow(gen);

        if( NewOrientation >= Orientation ){
            NewOrientation++;
        }

            L[s] = NewOrientation;
        }


        //Calcul de la contribution du spin modifié
        Matrix FacesFinal = Contact_Faces(L, NewLocation);
        for(int i=0;i<6;i++){
            if (FacesFinal(i,0) >= 0){
                DeltaE+=InteractionMap(FacesFinal(i,0), FacesFinal(i,1));
            }
        } 

        double ProbaAccept=exp(-beta*DeltaE);
    
        if (distribution(gen) < ProbaAccept){
            Particles(Particle_Number, 0) = NewLocation._x;
            Particles(Particle_Number, 1) = NewLocation._y;
        }
        else{
            L[NewLocation]=0;
            L[s]=Orientation;
            //On fait pas le changement
        }
}
int main() {
    
    int nx=10, ny=10, Nparts=2;
    float Taille=10;

    Lattice L=Lattice(nx,ny);
    Matrix Particles = Initialise_Lattice(L,Nparts);
    Matrix N_faces = initialise_face_count(L);

    Matrix InteractionMap = Matrix(6,6); 

    for(int i=0;i<6;i++){
        for(int j = 0; j<=i;j++){
            InteractionMap(i,j)= -10;
        }
    }

    //int steps = 0; //freqAffichage = pow(10,4);
    float beta=0.1;

    sf::RenderWindow window(sf::VideoMode(1000, 800), "Rotated Triangles inside Hexagon with SFML");
    sf::Color custom(127, 127, 127, 255); //Gris 

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear(custom); //Fond Gris
    
        Metropolis_step(L, Particles, InteractionMap, beta);


        L.affiche_SFML(window,Taille);
        
        window.display();
    }

}
