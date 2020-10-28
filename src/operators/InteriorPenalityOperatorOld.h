//
// Created by Alberto Chiappa on 21/12/17.
//

#ifndef PACS_PELI_CHIAPPA_INTERIORPENALITYOPERATOR_H
#define PACS_PELI_CHIAPPA_INTERIORPENALITYOPERATOR_H

#include <cmath>
#include "../problems/SystemMatrix.h"
#include "../domain/FeSpace.h"
#include "../basis/HierarchicalBasis.h"
#include "../basis/GradientHierarchicalBasis.h"
#include "../utilities/MatrixOperations.h"
#include "../utilities/GaussLobattoQuadrature.h"

class InteriorPenalityOperator {
private:
    FeSpace feSpace;
    Eigen::MatrixXd referenceIntegral;
    Eigen::MatrixXd referenceBoundaryLeft;
    Eigen::MatrixXd referenceBoundaryRight;
    Eigen::MatrixXd referenceBoundaryRightLeft;
    Eigen::MatrixXd referenceBoundaryLeftRight;


public:
    explicit InteriorPenalityOperator(FeSpace fespace);

    Eigen::MatrixXd computeReferenceIntegralMatrix();

    Eigen::MatrixXd computeReferenceRightBoundaryMatrix();

    Eigen::MatrixXd computeReferenceLeftBoundaryMatrix();

    Eigen::MatrixXd computeReferenceRightLeftBoundaryMatrix();

    Eigen::MatrixXd computeReferenceLeftRightBoundaryMatrix();

    void updateSystemMatrix(SystemMatrix &systemMatrix);

    int positionInCurrentElement(int i);

    int positionInAdjacentElement(int i);

    Eigen::MatrixXd computeLocalMatrix(const FeSpaceElement &firstElement, const FeSpaceElement &secondElement,
                                       int edgePosition, int otherEdgePosition);
};


#endif //PACS_PELI_CHIAPPA_INTERIORPENALITYOPERATOR_H
