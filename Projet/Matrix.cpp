#include"Matrix.h"

#include <iostream>
#include <stdexcept>
#include <fstream>

//_____________________________Fonctions

double Trace(const Matrix& M){
    double Tr=0;
    for(int k=0;k<std::min(M.nx,M.ny);k++){
        Tr+=M(k,k);
    }
    return Tr;
}

double Prodscal(const Matrix& A,const Matrix& B){
    Matrix At=A.T();
    return Trace(At*B);
}
//_____________________________Méthodes
void Matrix::afficher() const{
    
    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            std::cout << data[i*ny+j] << " ";
        };
        std::cout << "\n";
    };
}

Matrix& Matrix::T() const {
    Matrix* Mt = new Matrix(ny, nx);

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            (*Mt)(j, i) = (*this)(i, j);
        }
    }

    return *Mt;
}


//_____________________________Constructeur Matrix

Matrix::Matrix(int nx_, int ny_) : nx(nx_), ny(ny_), data(nullptr){

    data = new double[nx*ny];

    for(int k=0;k<nx*ny;k++){
        data[k]=0;
    }

}

Matrix::Matrix(const Matrix& other) : nx(other.nx), ny(other.ny), data(nullptr){
    data = new double[nx*ny];
    for (int k=0; k< nx*ny; k++){
        data[k] = other.data[k];
    }
}
//_____________________________Destructeur Matrix

Matrix::~Matrix(){
    if (data != nullptr){
        delete[] data;
    }
    data = nullptr;
}

//_____________________________Surcharge Operateurs Matrix

Matrix& Matrix::operator+(const Matrix& other) const{
    if (nx != other.nx || ny != other.ny) {
        throw std::invalid_argument("Matrices must have the same dimensions for addition.");
    }
    else{
        Matrix* sum= new Matrix(nx,ny);
        for(int k=0; k < nx*ny; k++){
            (*sum).data[k]=data[k]+other.data[k];
        }

        return *sum;
    }
}

Matrix& Matrix::operator-(const Matrix& other) const{
    if (nx != other.nx || ny != other.ny) {
        throw std::invalid_argument("Matrices must have the same dimensions for addition.");
    }
    else{
        Matrix* diff= new Matrix(nx,ny);
        for(int k=0; k < nx*ny; k++){
            (*diff).data[k]=data[k]-other.data[k];
        }

        return *diff;
    }
}

double Matrix::operator() (int i, int j) const{
    if((i > nx-1) || (i < 0) || (j > ny-1) || (j < 0) ){
        throw std::out_of_range("Matrix index out of range");
    }
    return data[i*ny+j];
    
}

double& Matrix::operator() (int i, int j){
    if((i > nx-1) || (i < 0) || (j > ny-1) || (j < 0) ){
        throw std::out_of_range("Matrix index out of range");
    }
    return data[i*ny+j];
}

Matrix& Matrix::operator+=(const Matrix& other){
    if (nx != other.nx || ny != other.ny) {
        throw std::invalid_argument("Matrices must have the same dimensions for addition.");
    }
    
    for(int k=0; k < nx*ny; k++){
        data[k]+=other.data[k];
    }

    return *this;
}

Matrix& Matrix::operator=(const Matrix& other){
    nx = other.nx;
    ny = other.ny;

    delete[] data;
    data = new double[nx*ny];

    for(int k=0; k < nx*ny; k++){
        data[k]=other.data[k];
    }

    return *this;
}

Matrix Matrix::operator*(const Matrix& other){
    if(ny != other.nx){
        throw std::invalid_argument("Matrices not compatible for multiplication.");
    }
    else{
        Matrix C = Matrix(nx,other.ny);

        for(int i=0;i<C.nx;i++){
            for(int j=0;j<C.ny;j++){
                for(int k=0;k<ny;k++){
                    C(i,j) += (*this)(i,k) * other(k,j);
                }
            }
        }

        return C;
    }
}

Matrix Matrix::operator*(const double& scalar){
    Matrix Result=Matrix(nx,ny);

    for(int i=0;i<nx;i++){
        for(int j=0; j<ny; j++){
            Result(i,j)= scalar * (*this)(i,j);
        }
    }

    return Result;
}

std::ostream& operator<<(std::ostream& os, const Matrix& M){
    os << "\n";
    for(int i=0; i< M.nx; i++){
        for(int j=0; j< M.ny; j++){
            std::cout << M(i,j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    return os;
}
//_____________________________Constructeur TriangularMatrix