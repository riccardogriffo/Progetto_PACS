#include "TransportEdgeOperator.h"

TransportEdgeOperator::TransportEdgeOperator(const FeSpace &fespace, double advectionCoeffX, double advectionCoeffY): fespace_(fespace),
                                                                                advectionCoeffX_(advectionCoeffX),advectionCoeffY_(advectionCoeffY) {
    unsigned int maxDegree = fespace.getMaxDegree();
    referenceElement_ = FeSpaceElement(maxDegree, maxDegree, std::make_pair(-1., 1.), std::make_pair(-1., 1.));
    computeReferenceMatrix1D();
}

void TransportEdgeOperator::computeReferenceMatrix1D() {
    unsigned int maxDegree = fespace_.getMaxDegree();
    referenceMatrix1D_ = Eigen::MatrixXd::Identity(maxDegree+1,maxDegree+1);
}

Eigen::MatrixXd TransportEdgeOperator::computeLocalMatrix(const FeSpaceElement &feSpaceElement, int edgeIndex,
                                                       const FeSpaceElement &adjacentElement) {

    int currentElemXDregree = feSpaceElement.getDegreeX();
    int currentElemYDregree = feSpaceElement.getDegreeY();
    int adjElemXDregree = adjacentElement.getDegreeX();
    int adjElemYDregree = adjacentElement.getDegreeY();

    Eigen::MatrixXd currentMatrix((currentElemXDregree + 1) * (currentElemYDregree + 1),
                                  (adjElemXDregree + 1) * (adjElemYDregree + 1));

    if (edgeIndex == 0) {

        // DA CAMBIARE
        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemYDregree + 1, feSpaceElement.getDomainY())
        );
        std::shared_ptr<HierarchicalBasis> basisPointerMinus(
                new HierarchicalBasis(adjElemYDregree + 1, feSpaceElement.getDomainY())
        );

        Eigen::VectorXd evalPoint(1);
        evalPoint << feSpaceElement.getDomainY().first;

        // valuto il jump della funzione test
        // auto evaluationPlus = basisPointerPlus->evaluate(evalPoint);
        // auto evaluationMinus = basisPointerMinus->evaluate(evalPoint);
        // auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);

        // calcolo il prodotto scalare nell'integrale
        // currentMatrix = MatrixOperations::tensorProduct(tempMatrix,
        //        referenceMatrix1D_.topLeftCorner(currentElemXDregree + 1,currentElemXDregree + 1));

        // Il calcolo del jump della funzione test rimane inalterato
        auto evaluationPlus = basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus = basisPointerMinus->evaluate(evalPoint);
        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);

        // nel calcolo del prodotto scalare devo sostiture il jump della soluzione
        // con il flusso numerico..
        // Qui si puo fare la versione che richiama un qualsiasi flusso numerico
        // Implemento ora Lax-Friedrics con beta costante (advectionCoeffX_, advectionCoeffY_)

        currentMatrix = -0.5*advectionCoeffX_ * MatrixOperations::tensorProduct(tempMatrix,1/feSpaceElement.getLengthX()*
                referenceMatrix1D_.topLeftCorner(currentElemXDregree + 1,currentElemXDregree + 1));

    }

    else if (edgeIndex == 2){

        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemYDregree+1, feSpaceElement.getDomainY())
        );

        std::shared_ptr<HierarchicalBasis> basisPointerMinus(
                new HierarchicalBasis(adjElemYDregree+1, adjacentElement.getDomainY())
        );
        Eigen::VectorXd evalPoint(1);
        evalPoint<<feSpaceElement.getDomainY().second;

        auto evaluationPlus=basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus=basisPointerMinus->evaluate(evalPoint);

        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);
        currentMatrix = 0.5*advectionCoeffX_ * MatrixOperations::tensorProduct(tempMatrix,1/feSpaceElement.getLengthX()*
                referenceMatrix1D_.topLeftCorner(currentElemXDregree + 1,currentElemXDregree + 1));
    }

    else if (edgeIndex == 1){

        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemXDregree + 1, feSpaceElement.getDomainX())
        );

        std::shared_ptr<HierarchicalBasis> basisPointerMinus(
                new HierarchicalBasis(adjElemXDregree + 1, adjacentElement.getDomainX())
        );
        Eigen::VectorXd evalPoint(1);
        evalPoint<<feSpaceElement.getDomainX().second;

        // valuto il jump della funzione test
        auto evaluationPlus=basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus=basisPointerMinus->evaluate(evalPoint);
        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);


        currentMatrix = -0.5*advectionCoeffY_ *MatrixOperations::tensorProduct(1/feSpaceElement.getLengthY()*
                                                  referenceMatrix1D_.topLeftCorner(currentElemYDregree + 1,
                                                                                   currentElemYDregree + 1),
                                                  tempMatrix);
    }

    else if (edgeIndex == 3){
        std::shared_ptr<HierarchicalBasis> basisPointerPlus(
                new HierarchicalBasis(currentElemXDregree + 1, feSpaceElement.getDomainX())
        );

        std::shared_ptr<HierarchicalBasis> basisPointerMinus(
                new HierarchicalBasis(adjElemXDregree + 1, adjacentElement.getDomainX())
        );

        Eigen::VectorXd evalPoint(1);
        evalPoint<<feSpaceElement.getDomainX().first;

        auto evaluationPlus=basisPointerPlus->evaluate(evalPoint);
        auto evaluationMinus=basisPointerMinus->evaluate(evalPoint);

        auto tempMatrix = MatrixOperations::tensorProduct(evaluationPlus.transpose(), evaluationMinus);
        currentMatrix = 0.5*advectionCoeffY_ *MatrixOperations::tensorProduct(1/feSpaceElement.getLengthY()*
                                                  referenceMatrix1D_.topLeftCorner(currentElemYDregree + 1,
                                                                                   currentElemYDregree + 1),
                                                  tempMatrix);
    }

    return currentMatrix;
}

void TransportEdgeOperator::updateSystemMatrix(SystemMatrix &systemMatrix) {

    auto fespaceElemVector = fespace_.getFeSpaceElements();
    for (int i = 0 ; i < fespaceElemVector.size(); ++i) {
        for (int j = 0; j < 4; ++j) {
            edge edge = fespace_.getMesh().getElement(i).getEdge(j);
            Eigen::MatrixXd currentMatrixUPlus = computeLocalMatrix(fespaceElemVector[i], j, fespaceElemVector[i]);
            if(edge.hasAdjacentElem()){
                Eigen::MatrixXd currentMatrixUMinus = computeLocalMatrix(fespaceElemVector[i], j,
                                                              fespaceElemVector[edge.getAdjacent()]);
                systemMatrix.addBlock(-currentMatrixUMinus, i, edge.getAdjacent());
                systemMatrix.addBlock(currentMatrixUPlus, i, i);
            }
            else{
                systemMatrix.addBlock(currentMatrixUPlus, i, i);
            }
        }
    }

}
