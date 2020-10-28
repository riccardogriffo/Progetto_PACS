//
// Created by Riccardo on 17/01/2018.
//

#include "Error.h"

Error::Error(std::function<double(double, double)> analiticalSolution, Eigen::VectorXd numericalSolution):
    analyticalSolution_(analiticalSolution),
    numericalSolution_(numericalSolution){
}

double Error::computeL2Error(FeSpace fespace) {
    RightHandSideFunction analyticalProjection (fespace, analyticalSolution_);
    Eigen::VectorXd analyticalCoeff = analyticalProjection.getSystemVector();

    Eigen::VectorXd coeffDiff = analyticalCoeff - numericalSolution_;

    return sqrt(coeffDiff.squaredNorm());
}

double Error::computeEnergyNorm(Eigen::VectorXd functionCoeff, const FeSpace &fespace, double gamma) {
    SystemMatrix systemMatrix(fespace);

    DiffusionOperator diffusionOperator(fespace, 1.);
    diffusionOperator.updateSystemMatrix(systemMatrix);

    StabilizerOperator stabilizerOperator(fespace, gamma);
    stabilizerOperator.updateSystemMatrix(systemMatrix);

    systemMatrix.buildSparseMatrix();

    auto matrix = systemMatrix.getMatrix();

    return sqrt(functionCoeff.adjoint()*matrix*functionCoeff);
}

double Error::computeEnergyError(FeSpace fespace, double gamma) {
    RightHandSideFunction analyticalProjection (fespace, analyticalSolution_);
    Eigen::VectorXd analyticalCoeff = analyticalProjection.getSystemVector();
    Eigen::VectorXd coeffDiff = analyticalCoeff - numericalSolution_;

    return computeEnergyNorm(coeffDiff, fespace, gamma);
}



