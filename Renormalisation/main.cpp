#include "TriangularLattice.h"
#include "Matrix.h"

#include <SFML/Graphics.hpp>
#include <cmath>
#include <iostream>
#include <algorithm> //min,max
#include <vector>
#include <fstream>

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

Matrix Pos_Particules(Lattice& L){
    int Nbr_particules= L.Particule_Count();    
    Matrix Particules = Matrix(Nbr_particules,2);

    int compteur = 0;
    for(int i=0;i<L.nx;i++){
        for(int j=0;j<L.ny;j++){
            if((int) L[L.site_xy(i,j)] != 0){
                Particules(compteur,0)=i;
                Particules(compteur,1)=j;
                compteur++;
            }
        }
    }

    return Particules;
}

Matrix Voisin_Face_Count(const Lattice& L, const Matrix& Particles){
    Matrix N_faces=Matrix(6,6);
    int site = 0, i, j;
    int face = 0;

    int min = 0;
    int max = 0;

    for(int k =0; k < Particles.nx;k++){
            i=Particles(k,0);
            j=Particles(k,1);
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
    return N_faces;
}

Matrix Second_Voisin_Face_Count(const Lattice& L, const Matrix& Particles){
    Matrix N_faces_second_voisin=Matrix(6,6);
    int site = 0, i, j;
    int face = 0;

    int min = 0;
    int max = 0;

    for(int k =0; k < Particles.nx;k++){
            i=Particles(k,0);
            j=Particles(k,1);
            site = (int) L[L.site_xy(i,j)];
            if(site != 0){
                face = 0;
                for(Site second_voisin: L.second_voisins(L.site_xy(i,j))){
                    if((int) L[second_voisin] != 0){
                        max = std::max((site-1+6-face)%6,( ((int) L[second_voisin])-1+3-face+6)%6); //On rajoute 6 pour que ça soit forcément positif
                        min = std::min((site-1+6-face)%6,( ((int) L[second_voisin])-1+3-face+6)%6);                       
                        N_faces_second_voisin(max,min) += 0.5;
                    }
                    face++;
                }
            }
    }
    return N_faces_second_voisin;
}

Matrix Mean_Face_Density(string Name, const int N_statistics){
    ifstream fich;
    fich.open(Name + ".dat");
    
    Matrix M=Matrix(6,6);
    double dens=0;

    for(int k=0;k<N_statistics;k++){
        for(int j=0;j<M.ny;j++){
            for(int i=j;i<M.nx;i++){
                fich >> dens;
                M(i,j)+=dens*(1./N_statistics);
            }
        }
    }

    fich.close();
    return M;
}

//Ecrit les entrées de la matrice triangulaire inférieure colonne par colonne
void Write_Matrix(const Matrix& M, ofstream& fich){
    for(int j=0;j<M.ny;j++){
        for(int i=j;i<M.nx;i++){
            fich << M(i,j) << " ";
        }
    }
    fich << "\n";
}

//__________________________Descente de Gradient
double eta(int t){
    return pow(10,2)*0.5*pow(0.97,t);
}

array<int,2> alpha_to_ij(int alpha){
    if((alpha < 0) || (alpha >20)){
        throw out_of_range("Alpha out of range");
    }
    else{
        return {alpha -5*(alpha>5) - 4*(alpha>10) - 3*(alpha>14) - 2*(alpha>17) - (alpha>19), (alpha>5) + (alpha>10) + (alpha>14) + (alpha>17) + (alpha>19)};
    }
}

int ij_to_alpha(int i,int j){
    if((j>i) || (i<0) || (i>5) || (j>5)){
        throw out_of_range("Indices out of range");
    }   
    else{
        return 5*j+i - (j>1) - 2*(j>2) - 3*(j>3) - 4*(j>4);
    }
}

double f(const double normc0squared, const Matrix& cm, const Matrix& d){
    return (1./(2*normc0squared))*Prodscal(cm-d,cm-d);
}

Matrix gradJf(const double normc0squared, const Matrix& cm, const Matrix& d, const Matrix& C){
    Matrix grad = Matrix(6,6);

    for(int beta=0;beta<C.ny;beta++){
        auto[ig,jg]=alpha_to_ij(beta);
        for(int alpha=0;alpha<C.nx;alpha++){
            auto[i,j]=alpha_to_ij(alpha);
            grad(ig,jg)+=(1/normc0squared)*(cm(i,j)-d(i,j))*C(max(alpha,beta), min(alpha,beta));
        }
    }
    return grad;
}

int main() {
    std::random_device rd;
    std::mt19937 gen(rd());

    int nx=256, ny=256;
    cin >> nx >> ny;
    int Nparts=(nx*ny)/9, N_bonds=3*nx*ny;
    //float Taille=8;

    Lattice L=Lattice(nx,ny);
    Matrix Particles = Initialise_Lattice(L,Nparts);
    //Matrix N_faces = Voisin_Face_Count(L, Particles);

    Matrix InteractionMap = Matrix(6,6); 
    
    //Matrix D=Matrix(6,6);
    ofstream fich2("Second_Voisins_Density_Faces.dat");

    //double mu = -2.5;
    //double sigma = 10.0;
    //std::normal_distribution<double> Normal(mu, sigma);
    //std::cout << Normal(gen);
    //for(int j=0;j<6;j++){
    //    for(int i=j;i<6;i++){
    //        InteractionMap(i,j) = Normal(gen);
    //    }
    //}

    InteractionMap(1,1)=-10, InteractionMap(3,2)=-10; InteractionMap(4,0)=-10;
    InteractionMap.afficher();
    cout << "___________\n";

    ofstream IntMap("Ancienne_Interaction_Map.txt");
    for(int i=0;i<InteractionMap.nx;i++){
        for(int j=0;j<InteractionMap.ny;j++){
            IntMap << InteractionMap(i,j) << " ";
        }
        IntMap << "\n";
    }
    IntMap.close();

    int N_T = 100; //Number of temperatures
    int N_steps = pow(10,4); //Number of steps by temperature
    //int N_annealing = N_T*N_steps;
    int N_statistics=pow(10,4); //Number of steps on which we do the averaging
    double beta = 1./N_T; //Température de départ

    int steps = 0;
    sf::Color custom(127, 127, 127, 255); //Gris 


    while(beta < 1){
        L.Metropolis_Step(L, Particles, InteractionMap, beta);   
        steps++;
        if ((steps%N_steps ==0) && (beta < 1)){
            beta+=1./N_T;
            steps=0;
        }
        else if (beta > 1){
            beta=1;
            std::cout << "kT=1 Atteint\n";
        }
    }
        
    for(int i=0;i<N_statistics;i++){
        L.Metropolis_Step(L, Particles, InteractionMap, 1);
        Write_Matrix(Second_Voisin_Face_Count(L, Particles)*(1./N_bonds),fich2);
        if(i%100==0){
        cout << "Statistique : " << i << "\n";
        }
    }
    fich2.close();


    Matrix d=Mean_Face_Density("Second_Voisins_Density_Faces", N_statistics);
    ofstream mat("Matrice_d.txt");
        for(int i=0;i<d.nx;i++){
            for(int j=0;j<d.ny;j++){
                mat << d(i,j) << " ";
            }
            mat << "\n";
        }
    mat.close();
    d.afficher();
    std::cout << "__________\n";

    Matrix NewInteractionMap = Matrix(6,6);
    Lattice R = L.Renormalisation();
    Matrix RParticles = Pos_Particules(R);
    int RN_bonds=3*R.nx*R.ny;

    //Matrix C=Matrix(6,6)
    ofstream fich1("Premier_Voisins_Density_Faces.dat");

    for(int k=0;k<N_steps;k++){
        R.Metropolis_Step(R, RParticles, NewInteractionMap, 1);
        Write_Matrix(Voisin_Face_Count(R, RParticles)*(1./RN_bonds), fich1);
    }
    fich1.close();

    Matrix c0=Mean_Face_Density("Premier_Voisins_Density_Faces", N_statistics);
    double normc0squared=Prodscal(c0,c0);
    Matrix v=Matrix(6,6);

    for(int t=0;t<100;t++){
        cout << "Decente de Gradient : " << t << "/100\n";
        Matrix cm = Mean_Face_Density("Premier_Voisins_Density_Faces", N_statistics);
        Matrix C = Matrix(21,21);
        Matrix c = Matrix(6,6);

        ifstream fich("Premier_Voisins_Density_Faces.dat");
        for(int k=0; k<N_statistics; k++){
            for(int j=0;j<c.ny;j++){
                for(int i=j;i<c.nx;i++){
                    fich >> c(i,j);
                }
            }
            for(int beta=0; beta<C.ny; beta++){
                for(int alpha=beta; alpha < C.nx; alpha++){
                    auto[i1,j1]=alpha_to_ij(alpha);
                    auto[i2,j2]=alpha_to_ij(beta);
                    C(alpha,beta)+= (c(i1,j1)-cm(i1,j1))*(c(i2,j2)-cm(i2,j2))*(1./N_statistics);
                }
            }
        }
        C=C*(-N_bonds);

        v = v*0.1 + gradJf(normc0squared, cm, d, C)*eta(t);
        NewInteractionMap = NewInteractionMap - v;

        ofstream newfich("Premier_Voisins_Density_Faces.dat");

        for(int k=0;k<N_steps;k++){
            R.Metropolis_Step(R, RParticles, NewInteractionMap, 1);
            Write_Matrix(Voisin_Face_Count(R, RParticles)*(1./RN_bonds), newfich);
        }
        if(t==99){
            c.afficher();
            ofstream mat("Matrice_c.txt");
            for(int i=0;i<c.nx;i++){
                for(int j=0;j<c.ny;j++){
                    mat << c(i,j) << " ";
                }
                mat << "\n";
            }

        }
            ofstream NewIntMap("Nouvelle_Interaction_Map.txt");
            for(int i=0;i<NewInteractionMap.nx;i++){
                for(int j=0;j<NewInteractionMap.ny;j++){
                    NewIntMap << NewInteractionMap(i,j) << " ";
                }
            NewIntMap << "\n";
            }
            NewIntMap.close();
    }
}
