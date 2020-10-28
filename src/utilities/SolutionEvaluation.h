//
// Created by Alberto Chiappa on 15/01/18.
//

#ifndef PACS_PELI_CHIAPPA_SOLUTIONEVALUATION_H
#define PACS_PELI_CHIAPPA_SOLUTIONEVALUATION_H


#include "../domain/FeSpace.h"
#include <Eigen/Dense>
#include <algorithm>

class SolutionEvaluation {
private:
    FeSpace fespace;

public:
    explicit SolutionEvaluation(FeSpace fespace);

    Eigen::MatrixXd evaluate(Eigen::VectorXd solutionCoeff, std::vector<double> xCoord, std::vector<double> yCoord);
};


#endif //PACS_PELI_CHIAPPA_SOLUTIONEVALUATION_H
