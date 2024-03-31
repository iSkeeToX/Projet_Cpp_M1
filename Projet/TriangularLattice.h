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

class Lattice{
    public:
    int nx,ny;

    void afficher() const;
    void affiche_SFML(sf::RenderWindow& window, float a) const;
    void write(std::string Name) const;

    Site site_xy(int x, int y) const;
    Site site_aleatoire() const;
    Site site_index(int index) const; 

    std::array<Site,6> voisins(const Site s) const; 
    std::array<Site,6> second_voisins(const Site s) const;
    
    int Particule_Count() const;
    Matrix Pos_Particules(const Lattice &L) const;

    int8_t operator[] (Site s) const;
    int8_t& operator[] (Site s);
    
    void Metropolis_Step(Lattice& L, Matrix& Particles, const Matrix& InteractionMap, const double beta);
    Lattice Renormalisation() const;
    
    //Crée un réseau triangulaire
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

    void write(const std::string Name) const;
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