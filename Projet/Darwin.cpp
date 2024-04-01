#include "Darwin.h"
#include "Matrix.h"

#include <iostream>
#include <random>


//_____________________________Le Constructeur de Darwin

Darwin::Darwin(int pop, int NbrGenes, double mean, double std) : pop(pop), NbrGenes(NbrGenes), Genes(Matrix(pop,NbrGenes)), Scores(Matrix(pop,2)){

    for(int i=0; i < Scores.nx; i++){
        Scores(i,1) = i;
    }

}

Darwin::~Darwin(){}