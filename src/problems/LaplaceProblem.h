//
// Created by Alberto Chiappa on 24/01/18.
//

#ifndef PACS_PELI_CHIAPPA_LAPLACEPROBLEM_H
#define PACS_PELI_CHIAPPA_LAPLACEPROBLEM_H


#include <iostream>
#include "../domain/FeSpace.h"
#include "../domain/QuadMesh.h"
#include "SystemMatrix.h"
#include "../operators/DiffusionOperator.h"
#include "../operators/InteriorPenalityOperator.h"
#include "../operators/StabilizerOperator.h"
#include "../utilities/SolutionEvaluation.h"
#include "../operators/RightHandSideFunction.h"
#include "../operators/InteriorPenalitySymmetric.h"
#include <Eigen/Dense>
#include <string>
#include "../utilities/Error.h"
#include <fstream>

class LaplaceProblem {
private:
    FeSpace fespace;
    double gamma;
    double tau;
    double diffCoeff;
    Eigen::VectorXd solutionCoeff;
    std::function<double (double, double)> f;
    bool solved;
public:
    LaplaceProblem(FeSpace fespace, double gamma, double tau, double diffCoeff,
                   std::function<double (double, double)> f);


    LaplaceProblem(std::string fileName, std::function<double (double, double)> f);
    void solve();

    void printOnScreen();

    void printOnFile(std::string fileName);

    double compareWithExactSolution(std::function<double (double, double)> exactSolution);

    double computeEnergyError(std::function<double (double, double)> exactSolution);
};


#endif //PACS_PELI_CHIAPPA_LAPLACEPROBLEM_H
