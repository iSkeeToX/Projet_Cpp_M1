#ifndef _TRIANG_LATTICE__
#define _TRIANG_LATTICE__

#include<cstdint>
#include <SFML/Graphics.hpp>

#include "Matrix.h"

class Site{
    public:
    int _index;
    int _x;
    int _y;

    friend class Lattice;
        public:

        Site (int i, int x, int y) : _index(i), _x(x), _y(y) {}
        Site () = default;
};

//La Classe Lattice est une classe nous permettant de créer un réseau triangulaire sur lequel siègent
//des hexagones réguliers qui peuvent se déplacer et s'orienter selon 6 directions (k\pi/3).
//
//Un site vide est repéré par un 0
//Un site occupé est repéré par un entier k compris entre 1 et 6 repérant son orientation (k\pi/3)
//par rapport à un axe horizontal orienté vers la droite
//
//Les données sont stockées dans un entier signé sur 8 bits permettant de réduire l'occupation de la mémoire
class Lattice{
    public:
    int nx,ny;

    void afficher() const;
    void affiche_SFML(sf::RenderWindow& window,const float a) const;

    Site site_xy(const int x,const int y) const;
    Site site_aleatoire() const;
    Site site_index(const int index) const; 

    std::array<Site,6> voisins(const Site s) const; 
    std::array<Site,6> second_voisins(const Site s) const;
    
    int Particule_Count() const;
    Matrix Pos_Particules(const Lattice &L) const;

    int8_t operator[] (Site s) const;
    int8_t& operator[] (Site s);
    
    void Metropolis_Step(Lattice& L, Matrix& Particles, const Matrix& InteractionMap, const double beta);
    Lattice Renormalisation() const;
    
    Lattice(int nx, int ny);
    ~Lattice();

    private:
    int8_t* data;
};

class SiteC{
    public:
    int _index;
    int _x;
    int _y;

    friend class ConComp;
        public:

        SiteC (int i, int x, int y) : _index(i), _x(x), _y(y) {}
        SiteC () = default;
};

class ConComp{
    public:
    int nx, ny, NbrCC;

    void Show_Connected_Components(const float a) const;

    SiteC site_xy(int x, int y) const;
    SiteC site_aleatoire() const;
    SiteC site_index(int index) const; 

    std::array<SiteC,6> voisins(const SiteC s) const; 

    int operator[] (SiteC s) const;
    int& operator[] (SiteC s);

    int OuterBorderLength(const int ConCompNumber) const;
    Matrix SizeConComps() const;
    
    bool IsOnLeftEdge(const int ConCompNumber) const;
    bool IsOnTopEdge(const int ConCompNumber) const;
    ConComp isolateConComp(const int ConCompNumber) const;

    ConComp Complementary() const;
    int SizeOfHoles() const;

    Matrix ClustersParameters() const;
    

    ConComp(Lattice&L);
    ConComp(int nx, int ny);
    ~ConComp();
    private:
    int* data;
};


#endif