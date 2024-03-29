#ifndef _MATRIX_
#define _MATRIX_

class Matrix{
    public:
    int nx,ny;

    void afficher() const;

    Matrix& T() const;

    Matrix(int nx,int ny);
    Matrix(const Matrix&);
    ~Matrix();

    double operator()(int i, int j) const;
    double& operator()(int i, int j);

    Matrix& operator+(const Matrix& other) const;
    Matrix& operator-(const Matrix& other) const;

    Matrix& operator+=(const Matrix& other);
    Matrix& operator=(const Matrix& other);
    
    Matrix operator*(const Matrix& other);
    Matrix operator*(const double& scalar);
    
    private:
    double* data;

};

class TriangularMatrix : public Matrix {





};

double Trace(const Matrix& M);
double Prodscal(const Matrix& A, const Matrix& B);


#endif