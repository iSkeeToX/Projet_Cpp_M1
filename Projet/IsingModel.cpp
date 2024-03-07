#include "IsingModel.h"

#include <random>
#include <iostream>
#include <algorithm> //min,max
#include <SFML/Graphics.hpp>
#include <SFML/System/Clock.hpp>

sf::Color custom(127, 127, 127, 255); //Gris 


//Initialise le lattice avec Nparticules ainsi que la matrice Nparts*2 contenant les coordonnées des particules (x,y)
void IsingModel::Initialise_Lattice(const int Nparts) {
    
    Particles = Matrix(Nparts,2);
    
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()  

    std::uniform_int_distribution<int8_t> orientation(1, 6);
    std::uniform_int_distribution<int> SiteAlea(0, nx*ny-1);
    
    int SiteAl=SiteAlea(gen);
    int8_t Or=orientation(gen);

    for(int i=0;i<Nparts;i++){
        SiteAl=SiteAlea(gen);
        if((*this)[(*this).site_index(SiteAl)]==0){
            Or=orientation(gen);
            (*this)[(*this).site_index(SiteAl)]=Or;

            Particles(i,0)= SiteAl/ny;//Coord x
            Particles(i,1)= (SiteAl + ny)%ny;
        }
        else{
            i--;
        }
    }

}

//Tire une Carte d'interaction suivant une loi normale
void IsingModel::Gaussian_InteractionMap(const float mean, const float standard_deviation){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::normal_distribution<double> distribution(mean, standard_deviation);

    for(int j=0;j<6;j++){
        for(int i=j; i<6;i++){
            InteractionMap(i,j) = distribution(gen);
        }
    }
}

//Met à jour la matrice contenant les positions des particules
void IsingModel::Pos_Particules(){

    int Nbr_particules= (*this).Particule_Count();    
    Particles = Matrix(Nbr_particules,2);

    int compteur = 0;
    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            if((int) (*this)[(*this).site_xy(i,j)] != 0){
                Particles(compteur,0)=i;
                Particles(compteur,1)=j;
                compteur++;
            }
        }
    }
}

//Renvoie la matrice 6*6 contenant le nombre d'occurence de chaque face (1er voisins)
Matrix IsingModel::Voisin_Face_Count() const{
    Matrix N_faces=Matrix(6,6);
    int site = 0, i, j;
    int face = 0;

    int minn = 0;
    int maxx = 0;

    for(int k =0; k < (*this).Particles.nx;k++){
            i=(*this).Particles(k,0);
            j=(*this).Particles(k,1);
            site = (int) (*this)[(*this).site_xy(i,j)];
            if(site != 0){
                face = 0;
                for(Site voisin: (*this).voisins((*this).site_xy(i,j))){
                    if((int) (*this)[voisin] != 0){
                        maxx = std::max((site-1+6-face)%6,( ((int) (*this)[voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
                        minn = std::min((site-1+6-face)%6,( ((int) (*this)[voisin])-1+3-face+6)%6);                       
                        N_faces(maxx,minn) += 0.5;
                    }
                    face++;
                }
            }
    }
    return N_faces;
}

//Renvoie la matrice 6*6 contenant le nombre d'occurence de chaque face (2nd voisins)
Matrix  IsingModel::Second_Voisin_Face_Count() const{
    Matrix N_faces_second_voisin=Matrix(6,6);
    int site = 0, i, j;
    int face = 0;

    int min = 0;
    int max = 0;

    for(int k =0; k < Particles.nx;k++){
            i=Particles(k,0);
            j=Particles(k,1);
            site = (int) (*this)[(*this).site_xy(i,j)];
            if(site != 0){
                face = 0;
                for(Site second_voisin: (*this).second_voisins((*this).site_xy(i,j))){
                    if((int) (*this)[second_voisin] != 0){
                        max = std::max((site-1+6-face)%6,( ((int) (*this)[second_voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
                        min = std::min((site-1+6-face)%6,( ((int) (*this)[second_voisin])-1+3-face+6)%6);                       
                        N_faces_second_voisin(max,min) += 0.5;
                    }
                    face++;
                }
            }
    }
    return N_faces_second_voisin;
}

Matrix IsingModel::Contact_Faces(const Site s, std::vector<int>& EmptySites) const{
    Matrix Faces = Matrix(6,2);
    int face=0;
    int max = 0,min = 0;

    int site = (int) (*this)[s];

    for(Site voisin: (*this).voisins(s)){
        if( (int) (*this)[voisin] == 0 ){
            Faces(face,0) = -1;
            Faces(face,1) = -1;
            EmptySites.push_back(face);
        }
        else{
            max = std::max((site-1+6-face)%6,( ((int) (*this)[voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
            min = std::min((site-1+6-face)%6,( ((int) (*this)[voisin])-1+3-face+6)%6);
            Faces(face,0)=max;
            Faces(face,1)=min;
        }
        face++;
    }
    return Faces;
}

Matrix IsingModel::Contact_Faces(const Site s) const{
    Matrix Faces = Matrix(6,2);
    int face=0;
    int max = 0,min = 0;

    int site = (int) (*this)[s];

    for(Site voisin: (*this).voisins(s)){
        if( (int) (*this)[voisin] == 0 ){
            Faces(face,0) = -1;
            Faces(face,1) = -1;
        }
        else{
            max = std::max((site-1+6-face)%6,( ((int) (*this)[voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
            min = std::min((site-1+6-face)%6,( ((int) (*this)[voisin])-1+3-face+6)%6);
            Faces(face,0)=max;
            Faces(face,1)=min;
        }
        face++;
    }
    return Faces;
}

void IsingModel::Metropolis_Step(){
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<int> RandomParticle(0, Particles.nx-1);
    
    int Particle_Number = RandomParticle(gen);
    Site s=(*this).site_xy(Particles(Particle_Number, 0), Particles(Particle_Number, 1));
    
    int Orientation = (*this)[s];
    double DeltaE=0;

    //Calcul de la contribution du site s à l'énergie
    std::vector<int> EmptySites; //Liste contenant l'indice des sites vides (ds l'ordre trig)  
    Matrix FacesInitial = (*this).Contact_Faces(s, EmptySites);
    
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
        std::uniform_int_distribution<int> MoveWhere(0,EmptySites.size()-1);
        NewLocation = (*this).voisins(s)[EmptySites[MoveWhere(gen)]];
        
        (*this)[s]=0;
        (*this)[NewLocation] = Orientation;      
    }
    else{
        std::uniform_int_distribution<int> OrientationHow(1,5);
        NewOrientation = OrientationHow(gen);

        if( NewOrientation >= Orientation ){
            NewOrientation++;
        }
        (*this)[s] = NewOrientation;
    }


    //Calcul de la contribution du spin modifié
    Matrix FacesFinal = (*this).Contact_Faces(NewLocation);
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
        (*this)[NewLocation]=0;
        (*this)[s]=Orientation;
        //On fait pas le changement
    }
}

void IsingModel::Annealing(const int N_T, const int N_steps, const float Taille, sf::RenderWindow& window){

    sf::Clock waitClock;
    beta = 1./N_T;

    int steps = 0;

    while((window.isOpen())){
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        
        window.clear(custom);
        for(int i = 0; (i < N_steps) && (steps < N_T); i++){
            (*this).Metropolis_Step();
        }
        (*this).affiche_SFML(window, Taille);
        window.display();
        
        if(steps < N_T){
        std::cout << "beta : " << beta << "\n______________________\n";
        beta+=1./N_T;
        steps++;
        }
    }
    

    beta = 1;
    std::cout << "Annealing terminé, beta = " << beta << "\n______________________\n";

}

void IsingModel::Annealing(const int N_T, const int N_steps){

    beta = 1./N_T;

    int steps = 0;
    while(steps < N_T){
        for(int i = 0; i < N_steps; i++){
            (*this).Metropolis_Step();
        }
        std::cout << "beta : " << beta << "\n______________________\n";
        steps++;
    }

    beta = 1;
    std::cout << "Annealing terminé, beta = " << beta << "\n______________________\n";

}

























