//
// Created by Alberto Chiappa on 15/01/18.
//

#include "SolutionEvaluation.h"

SolutionEvaluation::SolutionEvaluation(FeSpace fespace): fespace(fespace) {

}

Eigen::MatrixXd SolutionEvaluation::evaluate(Eigen::VectorXd solutionCoeff,
                                             std::vector<double> xCoord,
                                             std::vector<double> yCoord) {

    Eigen::MatrixXd solution = Eigen::MatrixXd::Zero(yCoord.size(), xCoord.size());
    auto elements = fespace.getFeSpaceElements();

    unsigned int elemIdx = 0;
    unsigned int currentSolutionIdx = 0;
    for (auto &elem : elements) {
        auto xDomain = elem.getDomainX();
        auto yDomain = elem.getDomainY();

        unsigned int xNumFunctions = elem.getDegreeX() + 1;
        unsigned int yNumFunctions = elem.getDegreeY() + 1;

        auto xFirstPtr = std::lower_bound(xCoord.begin(), xCoord.end(), xDomain.first);
        auto xSecondPtr = std::upper_bound(xCoord.begin(), xCoord.end(), xDomain.second);
        auto yFirstPtr = std::lower_bound(yCoord.begin(), yCoord.end(), yDomain.first);
        auto ySecondPtr = std::upper_bound(yCoord.begin(), yCoord.end(), yDomain.second);


        unsigned int numPointsX = xSecondPtr-xFirstPtr;
        unsigned int numPointsY = ySecondPtr-yFirstPtr;

        if (numPointsX > 0 && numPointsY >0) {
            Eigen::MatrixXd localSolution = Eigen::MatrixXd::Zero(numPointsY, numPointsX);
            Eigen::VectorXd localXCoord(numPointsX);
            Eigen::VectorXd localYCoord(numPointsY);

            for (auto xIt = xFirstPtr; xIt != xSecondPtr; ++xIt) {
                localXCoord[xIt-xFirstPtr] = *xIt;
            }
            for (auto yIt = yFirstPtr; yIt != ySecondPtr; ++yIt) {
                localYCoord[yIt-yFirstPtr] = *yIt;
            }

            Eigen::MatrixXd xBasisEval = elem.getBasisFunctionsX(localXCoord);
            Eigen::MatrixXd yBasisEval = elem.getBasisFunctionsY(localYCoord);

            for (int h = 0; h < numPointsX; ++h) {
                for (int k = 0; k < numPointsY; ++k) {
                    for (int i = 0; i <  xNumFunctions; ++i) {
                        for (int j = 0; j < yNumFunctions; ++j) {
                            localSolution(k, h) += solutionCoeff(currentSolutionIdx + i*yNumFunctions + j) *
                                                   xBasisEval(h, i) * yBasisEval(k, j);
                        }
                    }

                }
            }
            solution.block(yFirstPtr-yCoord.begin(), xFirstPtr-xCoord.begin(), numPointsY, numPointsX) = localSolution;
        }

        elemIdx++;
        currentSolutionIdx += xNumFunctions * yNumFunctions;
    }
    return solution;

}
