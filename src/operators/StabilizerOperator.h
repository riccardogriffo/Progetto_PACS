//
// Created by Riccardo on 17/12/2017.
//

#ifndef PACS_PELI_CHIAPPA_STABILIZEROPERATOR_H
#define PACS_PELI_CHIAPPA_STABILIZEROPERATOR_H

#include <Eigen/Dense>
#include "Operator.h"
#include "../domain/FeSpaceElement.h"
#include "../basis/HierarchicalBasis.h"
#include "IntegralMatrix1D.h"
#include "../utilities/MatrixOperations.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"


class StabilizerOperator: public Operator {
private:
    FeSpace fespace_;
    double gamma_;

public:
    StabilizerOperator(const FeSpace &fespace, double gamma);

    void computeReferenceMatrix1D() override;

    Eigen::MatrixXd computeLocalMatrix(const FeSpaceElement &feSpaceElement,
                                       int edgeIndex,
                                       const FeSpaceElement &adjacentElement);

    void updateSystemMatrix(SystemMatrix &systemMatrix);
};


#endif //PACS_PELI_CHIAPPA_STABILIZEROPERATOR_H
