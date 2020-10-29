#include "AdvectionOperator.h"

AdvectionOperator::AdvectionOperator(const FeSpace &fespace, double advectionCoeff): fespace_(fespace),
                                                                                     advectionCoeff_(advectionCoeff) {
    unsigned int maxDegree = fespace.getMaxDegree();
    referenceElement_ = FeSpaceElement(maxDegree, maxDegree, std::make_pair(-1., 1.), std::make_pair(-1., 1.));
    computeReferenceMatrix1D();
}

void AdvectionOperator::computeReferenceMatrix1D() {

    std::shared_ptr<GradientHierarchicalBasis> basisPointer_grad(
            new GradientHierarchicalBasis(referenceElement_.getDegreeX()+1, referenceElement_.getDomainX())
    );

    std::shared_ptr<HierarchicalBasis> basisPointer(
            new HierarchicalBasis(referenceElement_.getDegreeX()+1, referenceElement_.getDomainX())
    );

    IntegralMatrix1D matrix(basisPointer_grad,basisPointer);
    referenceMatrix1D_ = advectionCoeff_*matrix.generate(2*referenceElement_.getDegreeX());
}

Eigen::MatrixXd AdvectionOperator::computeLocalMatrix(FeSpaceElement feSpaceElement) {

    Eigen::MatrixXd identityX = Eigen::MatrixXd::Identity(feSpaceElement.getDegreeX()+1, feSpaceElement.getDegreeX()+1);
    Eigen::MatrixXd identityY = Eigen::MatrixXd::Identity(feSpaceElement.getDegreeY()+1, feSpaceElement.getDegreeY()+1);

    Eigen::MatrixXd referenceX = referenceMatrix1D_.topLeftCorner(feSpaceElement.getDegreeX()+1,feSpaceElement.getDegreeX()+1);
    Eigen::MatrixXd referenceY = referenceMatrix1D_.topLeftCorner(feSpaceElement.getDegreeY()+1,feSpaceElement.getDegreeY()+1);

    Eigen::MatrixXd firstAddend = MatrixOperations::tensorProduct(2/(feSpaceElement.getLengthX())*identityY,referenceX);
    Eigen::MatrixXd secondAddend = MatrixOperations::tensorProduct(2/(feSpaceElement.getLengthY())*referenceY,identityX);

    return firstAddend + secondAddend;
}

void AdvectionOperator::updateSystemMatrix(SystemMatrix &systemMatrix) {
    int curElemIdx = 0;
    auto vector = fespace_.getFeSpaceElements();
    for (auto &elem : vector) {
        Eigen::MatrixXd localAdvection = computeLocalMatrix(elem);
        systemMatrix.addBlock(localAdvection, curElemIdx, curElemIdx);
        ++curElemIdx;
