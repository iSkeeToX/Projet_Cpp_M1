#include<iostream>

#include "Matrix.h"


using namespace::std;

int main(){

    Matrix M = Matrix(2,2);

    M(0,1) = 1, M(1,0) = 2; M(1,1) = 3;

    cout << M;
}