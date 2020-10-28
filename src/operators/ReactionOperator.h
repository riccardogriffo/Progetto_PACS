//
// Created by Riccardo on 02/12/2017.
//

#ifndef PACS_PELI_CHIAPPA_REACTIONOPERATOR_H
#define PACS_PELI_CHIAPPA_REACTIONOPERATOR_H


#include "Operator.h"
#include "../utilities/MatrixOperations.h"
#include"../basis/HierarchicalBasis.h"
#include "IntegralMatrix1D.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"

class ReactionOperator: public Operator {
private:
    FeSpace fespace_;
    double reactionCoeff_;
public:
    ReactionOperator(const FeSpace &fespace, double reactionCoeff);

    void computeReferenceMatrix1D() override ;

    Eigen::MatrixXd computeLocalMatrix(FeSpaceElement feSpaceElement);

    void updateSystemMatrix(SystemMatrix &systemMatrix);

};


#endif //PACS_PELI_CHIAPPA_REACTIONOPERATOR_H
