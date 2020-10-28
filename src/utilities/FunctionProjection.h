//
// Created by Alberto Chiappa on 05/12/17.
//

#ifndef PACS_PELI_CHIAPPA_FUNCTIONPROJECTION_H
#define PACS_PELI_CHIAPPA_FUNCTIONPROJECTION_H


#include <Eigen/Dense>
#include "GaussLobattoQuadrature.h"
#include "../basis/Basis.h"
#include "MatrixOperations.h"
#include"../operators/IntegralMatrix1D.h"

#include <iostream>
#include <utility>
#include <memory>


class FunctionProjection {
private:
    GaussLobattoQuadrature quad;
    std::shared_ptr<Basis> xBasisPtr;
    std::shared_ptr<Basis> yBasisPtr;
    std::function<double (double, double)> func2d;

public:
    FunctionProjection(const std::function<double (double, double)> &func2d,
                       std::shared_ptr<Basis> xBasisPtr,
                       std::shared_ptr<Basis> yBasisPtr);

    Eigen::VectorXd computeProjection(unsigned int accuracyDegree);
};


#endif //PACS_PELI_CHIAPPA_FUNCTIONPROJECTION_H
