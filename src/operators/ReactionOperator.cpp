//
// Created by Riccardo on 02/12/2017.
//

#include "ReactionOperator.h"
#include "../domain/FeSpace.h"


ReactionOperator::ReactionOperator(const FeSpace &fespace, double reactionCoeff): fespace_(fespace),
                                                                                 reactionCoeff_(reactionCoeff) {
    unsigned int maxDegree = fespace.getMaxDegree();
    referenceElement_ = FeSpaceElement(maxDegree,maxDegree,std::make_pair(-1.,1.),std::make_pair(-1.,1.));
    computeReferenceMatrix1D();
}

void ReactionOperator::computeReferenceMatrix1D() {

    std::shared_ptr<HierarchicalBasis> basisPointer(
            new HierarchicalBasis(referenceElement_.getDegreeX()+1, referenceElement_.getDomainX())
    );

    IntegralMatrix1D matrix(basisPointer,basisPointer);
    referenceMatrix1D_ =reactionCoeff_*matrix.generate(2*(referenceElement_.getDegreeX()+1));
}

Eigen::MatrixXd ReactionOperator::computeLocalMatrix(FeSpaceElement feSpaceElement)  {
    return MatrixOperations::tensorProduct(referenceMatrix1D_.topLeftCorner(feSpaceElement.getDegreeX()+1,feSpaceElement.getDegreeX()+1),
                                     referenceMatrix1D_.topLeftCorner(feSpaceElement.getDegreeY()+1,feSpaceElement.getDegreeY()+1));
}

void ReactionOperator::updateSystemMatrix(SystemMatrix &systemMatrix) {
    int curElemIdx = 0;
    auto vector = fespace_.getFeSpaceElements();
    for (auto &elem : vector) {
        Eigen::MatrixXd localReaction = computeLocalMatrix(elem);
        systemMatrix.addBlock(localReaction, curElemIdx, curElemIdx);
        ++curElemIdx;
    }
}