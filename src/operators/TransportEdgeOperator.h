#ifndef PACS_PELI_CHIAPPA_TRANSPORTEDGEOPERATOR_H
#define PACS_PELI_CHIAPPA_TRANSPORTEDGEOPERATOR_H

#include <Eigen/Dense>
#include "Operator.h"
#include "../domain/FeSpaceElement.h"
#include "../basis/HierarchicalBasis.h"
#include "IntegralMatrix1D.h"
#include "../utilities/MatrixOperations.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"


class TransportEdgeOperator: public Operator {
private:
    FeSpace fespace_;
    double advectionCoeffX_;
    double advectionCoeffY_;


public:
    TransportEdgeOperator(const FeSpace &fespace, double advectionCoeffX, double advectionCoeffY);

    void computeReferenceMatrix1D() override;

    Eigen::MatrixXd computeLocalMatrix(const FeSpaceElement &feSpaceElement,
                                       int edgeIndex,
                                       const FeSpaceElement &adjacentElement);

    void updateSystemMatrix(SystemMatrix &systemMatrix);
};
