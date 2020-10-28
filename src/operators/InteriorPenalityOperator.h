//
// Created by Alberto Chiappa on 20/01/18.
//

#ifndef PACS_PELI_CHIAPPA_INTERIORPENALITYV2_H
#define PACS_PELI_CHIAPPA_INTERIORPENALITYV2_H


#include "Operator.h"
#include "../utilities/MatrixOperations.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"

class InteriorPenalityOperator: public Operator {
private:
    FeSpace fespace_;

public:
    InteriorPenalityOperator(const FeSpace &fespace);

    void computeReferenceMatrix1D() override;

    Eigen::MatrixXd computeLocalMatrix(const FeSpaceElement &feSpaceElement,
                                       int edgeIndex,
                                       const FeSpaceElement &adjacentElement);

    void updateSystemMatrix(SystemMatrix &systemMatrix);

};


#endif //PACS_PELI_CHIAPPA_INTERIORPENALITYV2_H
