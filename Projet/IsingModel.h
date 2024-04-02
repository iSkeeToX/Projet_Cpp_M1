#ifndef _ISING_MODEL_
#define _ISING_MODEL_

#include "TriangularLattice.h"

class IsingModel : public Lattice{
    public:
    Matrix InteractionMap;
    Matrix Particles;
    double beta;
    
    IsingModel(int nx, int ny) : Lattice(nx, ny), InteractionMap(Matrix(6,6)), Particles(Matrix(1,1)), beta(1) {}
    ~IsingModel() {};

    void Initialise_Lattice(const int Nparts);
    void Gaussian_InteractionMap(const float mean, const float standard_deviation);
    void Vider_Lattice();

    void Pos_Particules();
    Matrix Voisin_Face_Count() const;
    Matrix Second_Voisin_Face_Count() const;

    Matrix Contact_Faces(const Site s, std::vector<int>& EmptySites) const;
    Matrix Contact_Faces(const Site s) const;

    void Metropolis_Step();

    void Annealing(const int N_T, const int N_steps, const float Taille, sf::RenderWindow& window);
    void Annealing(const int N_T, const int N_steps);

    void TakePicture() const;
};


#endif