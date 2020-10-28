//
// Created by Alberto Chiappa on 24/11/17.
//

#ifndef PACS_PELI_CHIAPPA_OPERATORMATRIX_H
#define PACS_PELI_CHIAPPA_OPERATORMATRIX_H


#include "../basis/Basis.h"
#include "../utilities/GaussLobattoQuadrature.h"
#include <iostream>
#include <utility>
#include <memory>

class IntegralMatrix1D {

private:
    std::shared_ptr<Basis> firstBasisPtr;
    std::shared_ptr<Basis> secondBasisPtr;
    GaussLobattoQuadrature quad;
public:
    IntegralMatrix1D(std::shared_ptr<Basis> firstBasisPtr, std::shared_ptr<Basis> secondBasisPtr);
    Eigen::MatrixXd generate(unsigned int maxPolynomialDegree);

};


#endif //PACS_PELI_CHIAPPA_SYSTEMMATRIX_H
