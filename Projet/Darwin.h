#ifndef _DARWIN_H_
#define _DARWIN_H_

#include "Matrix.h"

class Darwin{
    public:
    
    int pop, generation, NbrGenes;
    Matrix Genes; //Genes(i,:) -> gènes du frérot i
    Matrix Scores; //Liste triée par fitness décroissante Genes(Scores(i,1),:) -> Gènes du ième meilleur

    Darwin(int pop, int NbrGenes, double mean, double std);
    ~Darwin();


};











#endif