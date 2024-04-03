#ifndef _MATRIX_
#define _MATRIX_

#include <iostream>

class Matrix{
    public:
    int nx,ny;

    void afficher() const;

    Matrix& T() const;
    
    void vstack(const Matrix& other);
    Matrix mean_columns() const;


    Matrix(int nx,int ny);
    Matrix(const Matrix&);
    Matrix(int nx, int ny, double mean, double std); 
    ~Matrix();

    double operator()(int i, int j) const;
    double& operator()(int i, int j);

    Matrix& operator+(const Matrix& other) const;
    Matrix& operator-(const Matrix& other) const;

    Matrix& operator+=(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    
    Matrix operator*(const Matrix& other);
    Matrix operator*(const double& scalar);
    
    friend std::ostream& operator<<(std::ostream& os, const Matrix& M);

    private:
    double* data;
};

double Trace(const Matrix& M);
double Prodscal(const Matrix& A, const Matrix& B);


#endif