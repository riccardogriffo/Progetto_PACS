//
// Created by Riccardo on 01/12/2017.
//

#include "DiffusionOperator.h"

DiffusionOperator::DiffusionOperator(const FeSpace &fespace, double diffusionCoeff): fespace_(fespace),
                                                                                     diffusionCoeff_(diffusionCoeff) {
    unsigned int maxDegree = fespace.getMaxDegree();
    referenceElement_ = FeSpaceElement(maxDegree, maxDegree, std::make_pair(-1., 1.), std::make_pair(-1., 1.));
    computeReferenceMatrix1D();
}

void DiffusionOperator::computeReferenceMatrix1D() {

    std::shared_ptr<GradientHierarchicalBasis> basisPointer(
            new GradientHierarchicalBasis(referenceElement_.getDegreeX()+1, referenceElement_.getDomainX())
    );
    IntegralMatrix1D matrix(basisPointer,basisPointer);
    referenceMatrix1D_ = diffusionCoeff_*matrix.generate(2*referenceElement_.getDegreeX());
}

Eigen::MatrixXd DiffusionOperator::computeLocalMatrix(FeSpaceElement feSpaceElement) {

    Eigen::MatrixXd identityX = Eigen::MatrixXd::Identity(feSpaceElement.getDegreeX()+1, feSpaceElement.getDegreeX()+1);
    Eigen::MatrixXd identityY = Eigen::MatrixXd::Identity(feSpaceElement.getDegreeY()+1, feSpaceElement.getDegreeY()+1);

    Eigen::MatrixXd referenceX = referenceMatrix1D_.topLeftCorner(feSpaceElement.getDegreeX()+1,feSpaceElement.getDegreeX()+1);
    Eigen::MatrixXd referenceY = referenceMatrix1D_.topLeftCorner(feSpaceElement.getDegreeY()+1,feSpaceElement.getDegreeY()+1);

    Eigen::MatrixXd firstAddend = MatrixOperations::tensorProduct(4/(feSpaceElement.getLengthX()*feSpaceElement.getLengthX())*
                                                                     identityY,referenceX);
    Eigen::MatrixXd secondAddend = MatrixOperations::tensorProduct(4/(feSpaceElement.getLengthY()*feSpaceElement.getLengthY())*
                                                                     referenceY,identityX);
    return firstAddend + secondAddend;
}

void DiffusionOperator::updateSystemMatrix(SystemMatrix &systemMatrix) {
    int curElemIdx = 0;
    auto vector = fespace_.getFeSpaceElements();
    for (auto &elem : vector) {
        Eigen::MatrixXd localDiffusion = computeLocalMatrix(elem);
        systemMatrix.addBlock(localDiffusion, curElemIdx, curElemIdx);
        ++curElemIdx;
    }
}




