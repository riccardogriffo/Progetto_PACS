

#ifndef PACS_PELI_CHIAPPA_DIFFUSIONOPERATOR_H
#define PACS_PELI_CHIAPPA_DIFFUSIONOPERATOR_H

#include "Operator.h"
#include "../utilities/MatrixOperations.h"
#include "../basis/GradientHierarchicalBasis.h"
#include "IntegralMatrix1D.h"
#include "../domain/FeSpace.h"
#include "../problems/SystemMatrix.h"

class AdvectionOperator: public Operator {
private:
    FeSpace fespace_;
    double advectionCoeff_;
public:

    AdvectionOperator(const FeSpace &fespace, double advectionCoeff);
    void computeReferenceMatrix1D() override;
    Eigen::MatrixXd computeLocalMatrix(FeSpaceElement feSpaceElement) ;
    void updateSystemMatrix(SystemMatrix &systemMatrix);

};


#endif //PACS_PELI_CHIAPPA_DIFFUSIONOPERATOR_H
