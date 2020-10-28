//
// Created by Riccardo on 17/01/2018.
//

#ifndef PACS_PELI_CHIAPPA_ERROR_H
#define PACS_PELI_CHIAPPA_ERROR_H

#include <Eigen/Dense>
#include <functional>
#include "../domain/FeSpace.h"
#include "../operators/RightHandSideFunction.h"
#include "../problems/SystemMatrix.h"
#include "../operators/DiffusionOperator.h"
#include "../operators/StabilizerOperator.h"


class Error {
private:
    std::function<double (double, double)> analyticalSolution_;
    Eigen::VectorXd numericalSolution_;
public:
    Error (std::function<double (double, double)> analyticalSolution, Eigen::VectorXd numericalSolution);
    double computeL2Error(FeSpace fespace);
    double computeEnergyError(FeSpace fespace, double gamma);
    double computeEnergyNorm(Eigen::VectorXd functionCoeff, const FeSpace &fespace, double gamma);

};


#endif //PACS_PELI_CHIAPPA_ERROR_H
