


#ifndef PACS_PELI_CHIAPPA_ADVECTIONDIFFUSIONROBLEM_H
#define PACS_PELI_CHIAPPA_ADVECTIONDIFFUSIONPROBLEM_H


#include <iostream>
#include "../domain/FeSpace.h"
#include "../domain/QuadMesh.h"
#include "SystemMatrix.h"
#include "../operators/AdvectionOperator.h"
#include "../operators/DiffusionOperator.h"
#include "../operators/InteriorPenalityOperator.h"
#include "../operators/TransportEdgeOperator.h"
#include "../operators/StabilizerOperator.h"
#include "../utilities/SolutionEvaluation.h"
#include "../operators/RightHandSideFunction.h"
#include "../operators/InteriorPenalitySymmetric.h"
#include <Eigen/Dense>
#include <string>
#include "../utilities/Error.h"
#include <fstream>

class AdvectionDiffusionProblem {
private:
    FeSpace fespace;
    double gamma;
    double tau;
    double diffCoeff;
    double advCoeffX;
    double advCoeffY;
    Eigen::VectorXd solutionCoeff;
    std::function<double (double, double)> f;
    bool solved;
public:
    AdvectionDiffusionProblem(FeSpace fespace, double gamma, double tau, double diffCoeff, double advCoeffX, double AdvCoeffY,
                   std::function<double (double, double)> f);


    AdvectionDiffusionProblem(std::string fileName, std::function<double (double, double)> f);
    void solve();

    void printOnScreen();

    void printOnFile(std::string fileName);

    double compareWithExactSolution(std::function<double (double, double)> exactSolution);

    double computeEnergyError(std::function<double (double, double)> exactSolution);
};


#endif //PACS_PELI_CHIAPPA_ADVECTIONDIFFUSIONROBLEM_H
