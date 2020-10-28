//
// Created by Alberto Chiappa on 12/12/17.
//

#include "RightHandSideFunction.h"
#include "../basis/HierarchicalBasis.h"
#include "../utilities/FunctionProjection.h"

RightHandSideFunction::RightHandSideFunction(FeSpace fespace,
                                             std::function<double(double, double)> func2d,
                                             unsigned int accuracyDegree): fespace(fespace),
                                                                           func2d(func2d),
                                                                           accuracyDegree(accuracyDegree) {
}

Eigen::VectorXd RightHandSideFunction::getSystemVector() {
    Eigen::VectorXd systemVector(fespace.getNumberBasisFunctions());
    std::vector<FeSpaceElement> elemVector = fespace.getFeSpaceElements();

    int sysVectorIdx = 0;

    for (auto elem = elemVector.begin(); elem != elemVector.end(); ++elem) {
        int xNumFunctions = elem->getDegreeX() + 1;
        int yNumFunctions = elem->getDegreeY() + 1;
        std::pair<double, double> xDomain = elem->getDomainX();
        std::pair<double, double> yDomain = elem->getDomainY();

        std::shared_ptr<Basis> xBasisPtr(new HierarchicalBasis(xNumFunctions, xDomain));
        std::shared_ptr<Basis> yBasisPtr(new HierarchicalBasis(yNumFunctions, yDomain));

        FunctionProjection functionProjection(func2d, xBasisPtr, yBasisPtr);

        Eigen::VectorXd localProjection = functionProjection.computeProjection(accuracyDegree);

        systemVector.segment(sysVectorIdx, xNumFunctions*yNumFunctions) = localProjection;

        sysVectorIdx += xNumFunctions*yNumFunctions;

    }
    return systemVector;
}
