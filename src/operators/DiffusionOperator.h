//
// Created by Riccardo on 01/12/2017.
//

#ifndef PACS_PELI_CHIAPPA_DIFFUSIONOPERATOR_H
#define PACS_PELI_CHIAPPA_DIFFUSIONOPERATOR_H

#include "Operator.h"
#include "../utilities/MatrixOperations.h"
#include "../basis/GradientHierarchicalBasis.h"
#include "IntegralMatrix1D.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"

class DiffusionOperator: public Operator {
private:
    FeSpace fespace_;
    double diffusionCoeff_;
public:

    DiffusionOperator(const FeSpace &fespace, double diffusionCoeff);
    void computeReferenceMatrix1D() override;
    Eigen::MatrixXd computeLocalMatrix(FeSpaceElement feSpaceElement) ;
    void updateSystemMatrix(SystemMatrix &systemMatrix);

};


#endif //PACS_PELI_CHIAPPA_DIFFUSIONOPERATOR_H
