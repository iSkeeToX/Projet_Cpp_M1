#include"TriangularLattice.h"

#include<iostream>
#include<random>
#include <stdexcept>
#include <algorithm> //min,max

std::random_device rd;  // Will be used to obtain a seed for the random number engine
std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

//_____________________________Fonctions nécessaires

Matrix Contact_Faces(const Lattice& L, const Site s, std::vector<int>& EmptySites){
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

//_____________________________Méthodes

void Lattice::afficher() const{
   for(int i=0; i<nx; i++){
    if(i%2==1){
        std::cout << "  ";
    }    
    for(int j=0; j<ny; j++){
        std::cout << (int) data[i*ny+j] << "   ";
    };
    std::cout << "\n";
    
   };
}

void hexagon(sf::RenderWindow& window, int8_t phi0, float a, float pos_x, float pos_y){

    sf::ConvexShape triangle;
    triangle.setPointCount(3);
    triangle.setPoint(0, sf::Vector2f(0, 0));
    triangle.setPoint(1, sf::Vector2f(a*cos(M_PI/6), a*sin(M_PI/6)));
    triangle.setPoint(2, sf::Vector2f(a*cos(M_PI/6), -a*sin(M_PI/6)));

    triangle.setPosition(pos_x, pos_y);
    triangle.setRotation((-phi0%6)*60);

    for(int i=0;i<6;i++){
        
        triangle.rotate(60);
        
        sf::Color triangleColor = sf::Color(30 * i, 50 * i, 70 * i);
        triangle.setFillColor(triangleColor);
        
        window.draw(triangle);
    }
}

void Lattice::affiche_SFML(sf::RenderWindow& window, float a) const{
    float taille=a*sqrt(3)/2;

    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            int8_t valeur=data[i*ny+j];
            if(valeur!=0){
                hexagon(window,valeur,a,taille*(2*j+1+i%2),a*(1.5*i+1));
            }
        }        
    }
}

Site Lattice::site_xy(int x, int y) const{
    int xm=x%nx;
    int ym=y%ny;
    if (xm < 0 && ym < 0){
        xm+=nx; ym+=ny;
        return Site(ny*xm+ym,xm,ym);
    }
    else if (xm < 0){
        xm+=nx;
        return Site(ny*xm+ym,xm,ym);
    }
    else if (ym < 0){
        ym+=ny;
        return Site(ny*xm+ym,xm,ym);
    }
    else{
        return Site(ny*xm+ym,xm,ym);
    }
}

Site Lattice::site_aleatoire() const{
    
    std::uniform_int_distribution<int> distribution(0, nx*ny);
    
    int _index=distribution(gen);
    int y=_index%ny;
    if(y < 0){
        return Site(_index,_index/ny,y+ny);
    }
    else{
        return Site(_index,_index/ny,y);
    }
}

//index=ny*x+y
Site Lattice::site_index(int index) const{
    int y=index%ny;
    if(y < 0){
        return Site(index,index/ny,y+ny);
    }
    else{
        return Site(index,index/ny,y);
    }
}

//Tourne dans le sens trigo, en partant de la droite
std::array<Site,6> Lattice::voisins(const Site s) const{
    int x=s._x, y=s._y;
    if(x%2 ==0){
        return {this->site_xy(x,y+1), this->site_xy(x-1,y), this->site_xy(x-1,y-1), this->site_xy(x,y-1), this->site_xy(x+1,y-1), this->site_xy(x+1,y)};
    }
    else{
        return {this->site_xy(x,y+1), this->site_xy(x-1,y+1), this->site_xy(x-1,y), this->site_xy(x,y-1), this->site_xy(x+1,y), this->site_xy(x+1,y+1)};
    }
}

//Tourne dans le sens trigo, en partant de la droite
std::array<Site,6> Lattice::second_voisins(const Site s) const{
    int x=s._x, y=s._y;
    if(x%2 ==0){
        return {this->site_xy(x,y+2), this->site_xy(x-2,y+1), this->site_xy(x-2,y-1), this->site_xy(x,y-2), this->site_xy(x+2,y-1), this->site_xy(x+2,y+1)};
    }
    else{
        return {this->site_xy(x,y+2), this->site_xy(x-2,y+1), this->site_xy(x-2,y-1), this->site_xy(x,y-2), this->site_xy(x+2,y-1), this->site_xy(x+2,y+1)};
    }
}

int Lattice::Particule_Count() const{
    int count =0;

    for(int k=0;k<nx*ny;k++){
        if(data[k]!=0){
            count+=1;
        }
    }
    return count;
}

void Lattice::Metropolis_Step(Lattice& L, Matrix& Particles, const Matrix& InteractionMap, const double beta){

        std::uniform_int_distribution<int> RandomParticle(0, Particles.nx-1);
        int Particle_Number = RandomParticle(gen);
        Site s=L.site_xy(Particles(Particle_Number, 0), Particles(Particle_Number, 1));
    
        int Orientation = L[s];
        double DeltaE=0;

        //Calcul de la contribution du site s à l'énergie
        std::vector<int> EmptySites; //Liste contenant l'indice des sites vides (ds l'ordre trig)  
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
            std::uniform_int_distribution<int> MoveWhere(0,EmptySites.size()-1);
            NewLocation = L.voisins(s)[EmptySites[MoveWhere(gen)]];
        
            L[s]=0;
            L[NewLocation] = Orientation;      
        }
        else{
            std::uniform_int_distribution<int> OrientationHow(1,5);
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

Lattice Lattice::Renormalisation() const{
    Lattice R = Lattice(nx/2,ny/2);

    for(int i=0;i<nx;i++){
        for(int j=0;j<ny;j++){
            if((i%2==0) && (j%2==0)){
                if(i%4==2){
                    R[R.site_xy(i/2,j/2)] = (*this)[this->site_xy(i,j+1)];
                }
                else{
                    R[R.site_xy(i/2,j/2)] = (*this)[this->site_xy(i,j)];
                }
            }
        }
    }
    return R;
}
//_____________________________Opérateurs

int8_t Lattice::operator[] (Site s) const{
    return data[s._index];
}

int8_t& Lattice::operator[] (Site s){
    return data[s._index];
}


//_____________________________Constructeurs

Lattice::Lattice(int nx_, int ny_) : nx(nx_), ny(ny_), data(nullptr){
    data= new int8_t[nx*ny];
    
    if(ny%2==1){
        throw std::runtime_error("ny doit être pair pour respecter les conditions de BVK");
    }
    for(int i=0;i<nx*ny;i++){
        data[i]=0;
    }
}

//Destructeur
Lattice::~Lattice(){
    if (data != nullptr){
        delete[] data;
    }
    data = nullptr;
}

//_____________________________Privée

