//
// Created by Riccardo on 21/01/2018.
//

#ifndef PACS_PELI_CHIAPPA_INTERIORPENALITYSYMMETRIC_H
#define PACS_PELI_CHIAPPA_INTERIORPENALITYSYMMETRIC_H

#include "Operator.h"
#include "../utilities/MatrixOperations.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"

class InteriorPenalitySymmetric : public Operator {
private:
    FeSpace fespace_;
    double tau_;

public:
    InteriorPenalitySymmetric(const FeSpace &fespace, double tau);

    void computeReferenceMatrix1D() override;

    Eigen::MatrixXd computeLocalMatrix(const FeSpaceElement &feSpaceElement,
                                       int edgeIndex,
                                       const FeSpaceElement &adjacentElement);

    void updateSystemMatrix(SystemMatrix &systemMatrix);

};


#endif //PACS_PELI_CHIAPPA_INTERIORPENALITYSYMMETRIC_H
