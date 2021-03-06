//
// Created by Alberto Chiappa on 20/01/18.
//

#include "InteriorPenalityOperator.h"

InteriorPenalityOperator::InteriorPenalityOperator(const FeSpace &fespace): fespace_(fespace) {
    unsigned int maxDegree = fespace.getMaxDegree();
    referenceElement_ = FeSpaceElement(maxDegree, maxDegree, std::make_pair(-1., 1.), std::make_pair(-1., 1.));
    computeReferenceMatrix1D();
}

void InteriorPenalityOperator::computeReferenceMatrix1D() {
    unsigned int maxDegree = fespace_.getMaxDegree();
    referenceMatrix1D_ =Eigen::MatrixXd::Identity(maxDegree+1,maxDegree+1);
}

Eigen::MatrixXd InteriorPenalityOperator::computeLocalMatrix(const FeSpaceElement &feSpaceElement, int edgeIndex,
                                                       const FeSpaceElement &adjacentElement) {
    int currentElemXDregree = feSpaceElement.getDegreeX();
    int currentElemYDregree = feSpaceElement.getDegreeY();
    int adjElemXDregree = adjacentElement.getDegreeX();
    int adjElemYDregree = adjacentElement.getDegreeY();

    Eigen::MatrixXd currentMatrix((currentElemXDregree + 1) * (currentElemYDregree + 1),
                                  (adjElemXDregree + 1) * (adjElemYDregree + 1));

    if (edgeIndex == 0) {

        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemYDregree + 1, feSpaceElement.getDomainY())
        );

        std::shared_ptr<GradientHierarchicalBasis> basisPointerMinus(
                new GradientHierarchicalBasis(adjElemYDregree + 1, adjacentElement.getDomainY())
        );
        Eigen::VectorXd evalPoint(1);
        evalPoint << feSpaceElement.getDomainY().first;

        auto evaluationPlus = basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus = basisPointerMinus->evaluate(evalPoint);

        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);

        currentMatrix = -0.5*MatrixOperations::tensorProduct(tempMatrix,referenceMatrix1D_.topLeftCorner(currentElemXDregree + 1,currentElemXDregree + 1));
    }

    else if (edgeIndex == 2){

        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemYDregree+1, feSpaceElement.getDomainY())
        );

        std::shared_ptr<GradientHierarchicalBasis> basisPointerMinus(
                new GradientHierarchicalBasis(adjElemYDregree+1, adjacentElement.getDomainY())
        );
        Eigen::VectorXd evalPoint(1);
        evalPoint<<feSpaceElement.getDomainY().second;

        auto evaluationPlus=basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus=basisPointerMinus->evaluate(evalPoint);

        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);
        currentMatrix = 0.5*MatrixOperations::tensorProduct(tempMatrix,referenceMatrix1D_.topLeftCorner(currentElemXDregree + 1,currentElemXDregree + 1));
    }

    else if (edgeIndex == 1){

        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemXDregree + 1, feSpaceElement.getDomainX())
        );

        std::shared_ptr<GradientHierarchicalBasis> basisPointerMinus(
                new GradientHierarchicalBasis(adjElemXDregree + 1, adjacentElement.getDomainX())
        );
        Eigen::VectorXd evalPoint(1);
        evalPoint<<feSpaceElement.getDomainX().second;

        auto evaluationPlus=basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus=basisPointerMinus->evaluate(evalPoint);

        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);
        currentMatrix = 0.5*MatrixOperations::tensorProduct(referenceMatrix1D_.topLeftCorner(currentElemYDregree + 1,
                                                                                   currentElemYDregree + 1),
                                                  tempMatrix);
    }

    else if (edgeIndex == 3){
        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemXDregree + 1, feSpaceElement.getDomainX())
        );

        std::shared_ptr<GradientHierarchicalBasis> basisPointerMinus(
                new GradientHierarchicalBasis(adjElemXDregree + 1, adjacentElement.getDomainX())
        );

        Eigen::VectorXd evalPoint(1);
        evalPoint<<feSpaceElement.getDomainX().first;

        auto evaluationPlus=basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus=basisPointerMinus->evaluate(evalPoint);

        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);
        currentMatrix = -0.5*MatrixOperations::tensorProduct(referenceMatrix1D_.topLeftCorner(currentElemYDregree + 1,
                                                                                   currentElemYDregree + 1),
                                                  tempMatrix);
    }

    return currentMatrix;
}

void InteriorPenalityOperator::updateSystemMatrix(SystemMatrix &systemMatrix) {
    auto fespaceElemVector = fespace_.getFeSpaceElements();
    for (int i = 0 ; i < fespaceElemVector.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            auto edge = fespace_.getMesh().getElement(i).getEdge(j);
            auto currentMatrixUPlus = computeLocalMatrix(fespaceElemVector[i], j, fespaceElemVector[i]);
            systemMatrix.addBlock(-currentMatrixUPlus, i, i);
            if(edge.hasAdjacentElem()){
                auto currentMatrixUMinus = computeLocalMatrix(fespaceElemVector[i], j,
                                                              fespaceElemVector[edge.getAdjacent()]);
                systemMatrix.addBlock(-currentMatrixUMinus, i, edge.getAdjacent());
            }
            else{
                systemMatrix.addBlock(-currentMatrixUPlus, i, i);
            }
        }
    }
}
