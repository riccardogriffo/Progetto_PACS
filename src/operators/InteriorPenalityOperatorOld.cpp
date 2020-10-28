//
// Created by Alberto Chiappa on 21/12/17.
//

#include "InteriorPenalityOperatorOld.h"

#define BOTTOM 0
#define TOP 2
#define LEFT 1
#define RIGHT 3

InteriorPenalityOperator::InteriorPenalityOperator(FeSpace fespace): feSpace(fespace) {
    this->referenceIntegral = this->computeReferenceIntegralMatrix();
    this->referenceBoundaryLeft = this->computeReferenceLeftBoundaryMatrix();
    this->referenceBoundaryRight = this->computeReferenceRightBoundaryMatrix();
    this->referenceBoundaryRightLeft = this->computeReferenceRightLeftBoundaryMatrix();
    this->referenceBoundaryLeftRight = this->computeReferenceLeftRightBoundaryMatrix();

}

void InteriorPenalityOperator::updateSystemMatrix(SystemMatrix &systemMatrix) {

    quadMesh mesh = this->feSpace.getMesh();
    std::vector<quadElement> elements = mesh.getElements();
    std::vector<FeSpaceElement> feSpaceElementVec = feSpace.getFeSpaceElements();
    Eigen::MatrixXd localPenality;

    int curElemIdx = 0;
    for (auto &element : elements) {
        FeSpaceElement feSpaceElement = feSpaceElementVec[curElemIdx];

        std::vector<edge> edges = element.getEdges();
        for (int i = 0; i < 4; ++i) {
            localPenality = computeLocalMatrix(feSpaceElement, feSpaceElement, positionInCurrentElement(i),
                                               positionInCurrentElement(i));
            systemMatrix.addBlock(-localPenality, curElemIdx, curElemIdx);

            edge edge = edges[i];
            if (edge.hasAdjacentElem()) {
                int otherElemIdx = edge.getAdjacent();
                FeSpaceElement otherElement = feSpaceElementVec[otherElemIdx];
                localPenality = computeLocalMatrix(feSpaceElement, otherElement, positionInCurrentElement(i),
                                                   positionInAdjacentElement(i));
                systemMatrix.addBlock(-localPenality, curElemIdx, otherElemIdx);
            }
            else {
                systemMatrix.addBlock(-localPenality, curElemIdx, curElemIdx);
            }
        }
        ++curElemIdx;
    }
}

Eigen::MatrixXd InteriorPenalityOperator::computeReferenceIntegralMatrix() {
    GaussLobattoQuadrature quad;

    unsigned int maxBasisDegree = feSpace.getMaxDegree();
    unsigned int maxIntegrationDegree = 2*maxBasisDegree;
    unsigned int numNodes = (maxIntegrationDegree + 3)/2 + (maxIntegrationDegree % 2 == 0);

    HierarchicalBasis hierarchicalBasis(maxBasisDegree + 1);

    Eigen::VectorXd nodes;
    Eigen::VectorXd weights;
    std::tie(nodes, weights) = quad.computeNodesAndWeights(numNodes, hierarchicalBasis.getDomain());

    Eigen::MatrixXd basisEval = hierarchicalBasis.evaluate(nodes);

    Eigen::MatrixXd referenceMatrix(maxBasisDegree + 1, maxBasisDegree + 1);

    for (int i = 0; i < maxBasisDegree + 1; ++i) {
        for (int j = 0; j < maxBasisDegree + 1; ++j) {
            double result = MatrixOperations().vectorElementWiseProduct(basisEval.col(i), basisEval.col(j)).dot(weights);
            referenceMatrix(i, j) = result;
        }
    }

    return referenceMatrix;
}

/*
 *  The index of the row corresponds to the degree of the basis polynomial, the index of column corresponds to the
 *  degree of the gradient of the basis polynomial -1.
 */
Eigen::MatrixXd InteriorPenalityOperator::computeReferenceRightBoundaryMatrix() {
    Eigen::VectorXd boundary(1);
    boundary(0) = 1.;

    unsigned int maxBasisDegree = feSpace.getMaxDegree();

    HierarchicalBasis hierarchicalBasis(maxBasisDegree + 1) ;
    Eigen::MatrixXd boundaryEval = hierarchicalBasis.evaluate(boundary);

    GradientHierarchicalBasis gradientHierarchicalBasis(maxBasisDegree + 1);
    Eigen::MatrixXd gradBoundaryEval = gradientHierarchicalBasis.evaluate(boundary);

    Eigen::MatrixXd referenceMatrix = MatrixOperations().tensorProduct(boundaryEval.transpose(), gradBoundaryEval);

    return referenceMatrix;
}

/*
 *  The index of the row corresponds to the degree of the basis polynomial, the index of column corresponds to the
 *  degree of the gradient of the basis polynomial -1.
 */
Eigen::MatrixXd InteriorPenalityOperator::computeReferenceLeftBoundaryMatrix() {
    Eigen::VectorXd boundary(1);
    boundary(0) = -1.;

    unsigned int maxBasisDegree = feSpace.getMaxDegree();

    HierarchicalBasis hierarchicalBasis(maxBasisDegree + 1) ;
    Eigen::MatrixXd boundaryEval = hierarchicalBasis.evaluate(boundary);

    GradientHierarchicalBasis gradientHierarchicalBasis(maxBasisDegree + 1);
    Eigen::MatrixXd gradBoundaryEval = gradientHierarchicalBasis.evaluate(boundary);

    Eigen::MatrixXd referenceMatrix = MatrixOperations().tensorProduct(boundaryEval.transpose(), gradBoundaryEval);

    return referenceMatrix;
}

Eigen::MatrixXd InteriorPenalityOperator::computeReferenceRightLeftBoundaryMatrix() {
    Eigen::VectorXd boundaryRight(1);
    boundaryRight(0) = 1.;

    Eigen::VectorXd boundaryLeft(1);
    boundaryLeft(0) = -1.;

    unsigned int maxBasisDegree = feSpace.getMaxDegree();

    HierarchicalBasis hierarchicalBasis(maxBasisDegree + 1) ;
    Eigen::MatrixXd boundaryEval = hierarchicalBasis.evaluate(boundaryRight);

    GradientHierarchicalBasis gradientHierarchicalBasis(maxBasisDegree + 1);
    Eigen::MatrixXd gradBoundaryEval = gradientHierarchicalBasis.evaluate(boundaryLeft);

    Eigen::MatrixXd referenceMatrix = MatrixOperations().tensorProduct(boundaryEval.transpose(), gradBoundaryEval);

    return referenceMatrix;
}

Eigen::MatrixXd InteriorPenalityOperator::computeReferenceLeftRightBoundaryMatrix() {
    Eigen::VectorXd boundaryRight(1);
    boundaryRight(0) = 1.;

    Eigen::VectorXd boundaryLeft(1);
    boundaryLeft(0) = -1.;

    unsigned int maxBasisDegree = feSpace.getMaxDegree();

    HierarchicalBasis hierarchicalBasis(maxBasisDegree + 1) ;
    Eigen::MatrixXd boundaryEval = hierarchicalBasis.evaluate(boundaryLeft);

    GradientHierarchicalBasis gradientHierarchicalBasis(maxBasisDegree + 1);
    Eigen::MatrixXd gradBoundaryEval = gradientHierarchicalBasis.evaluate(boundaryRight);

    Eigen::MatrixXd referenceMatrix = MatrixOperations().tensorProduct(boundaryEval.transpose(), gradBoundaryEval);

    return referenceMatrix;}

int InteriorPenalityOperator::positionInCurrentElement(int edgeIdx) {
    switch (edgeIdx) {
        case 0:
            return BOTTOM;
        break;

        case 1:
            return RIGHT;
        break;

        case 2:
            return TOP;
        break;

        case 3:
            return LEFT;
        break;

        default:
            throw std::invalid_argument("The accepted value are only 0, 1, 2 or 3.");
        break;
    }
}

int InteriorPenalityOperator::positionInAdjacentElement(int edgeIdx) {
    switch (edgeIdx) {
        case 0:
            return TOP;
        break;

        case 1:
            return LEFT;
        break;

        case 2:
            return BOTTOM;
        break;

        case 3:
            return RIGHT;
        break;

        default:
            throw std::invalid_argument("The accepted value are only 0, 1, 2 or 3.");
        break;
    }
}

/*
 * The first element is the one that states on which column of the big matrix we are writing, because it states which
 * coefficient of u has to be multuplied. The second element identifies the basis function v, therfore statin on which
 * column we are going to write. The first element indices also pick the values of the derivative of the basis function
 * at the boundary, while the indeces of the second one pick the non derived values.
 */
Eigen::MatrixXd
InteriorPenalityOperator::computeLocalMatrix(const FeSpaceElement &firstElement, const FeSpaceElement &secondElement, int edgePosition,
                                     int otherEdgePosition) {
    unsigned int xDegreeFirst = firstElement.getDegreeX();
    unsigned int yDegreeFirst = firstElement.getDegreeY();
    unsigned int xDegreeSecond = secondElement.getDegreeX();
    unsigned int yDegreeSecond = secondElement.getDegreeY();

    double xDomain = secondElement.getLengthX();
    double yDomain = secondElement.getLengthY();

    double xDomainFirst;
    double yDomainFirst;
    double xDomainSecond;
    double yDomainSecond;

    Eigen::MatrixXd localPenality((xDegreeFirst+1)*(yDegreeFirst+1),(xDegreeSecond+1)*(yDegreeSecond+1));
    if (edgePosition == otherEdgePosition) {
        switch (edgePosition) {
            case RIGHT:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = 0.5 * referenceBoundaryRight(h, i) * referenceIntegral(k, j);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                xDomainFirst = firstElement.getLengthX();
                localPenality *= 4/xDomainFirst/xDomainFirst;
            break;

            case TOP:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = 0.5 * referenceBoundaryRight(k, j) * referenceIntegral(h, i);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                yDomainFirst = firstElement.getLengthY();
                localPenality *= 4/yDomainFirst/yDomainFirst;
            break;

            case LEFT:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = -0.5 * referenceBoundaryLeft(h, i) * referenceIntegral(k, j);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                xDomainFirst = firstElement.getLengthX();
                localPenality *= 4/xDomainFirst/xDomainFirst;
                break;

            case BOTTOM:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = -0.5 * referenceBoundaryLeft(k, j) * referenceIntegral(h, i);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                yDomainFirst = firstElement.getLengthY();
                localPenality *= 4/yDomainFirst/yDomainFirst;
            break;

            default:
                throw std::invalid_argument("The accepted value are only 0, 1, 2 or 3.");
                break;
        }
    }
    else {
        switch (otherEdgePosition) {
            case RIGHT:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = -0.5 * referenceBoundaryLeftRight(h, i) * referenceIntegral(k, j);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                xDomainFirst = firstElement.getLengthX();
                xDomainSecond = secondElement.getLengthX();
                localPenality *= 4/sqrt(xDomainFirst)/sqrt(xDomainSecond)/xDomainSecond;
            break;

            case TOP:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = -0.5 * referenceBoundaryLeftRight(k, j) * referenceIntegral(h, i);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                yDomainFirst = firstElement.getLengthY();
                yDomainSecond = secondElement.getLengthY();
                localPenality *= 4/sqrt(yDomainFirst)/sqrt(yDomainSecond)/yDomainSecond;
            break;

            case LEFT:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = 0.5 * referenceBoundaryRightLeft(h, i) * referenceIntegral(k, j);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                xDomainFirst = firstElement.getLengthX();
                xDomainSecond = secondElement.getLengthX();
                localPenality *= 4/sqrt(xDomainFirst)/sqrt(xDomainSecond)/xDomainSecond;
            break;

            case BOTTOM:
                for (int i = 0; i < xDegreeSecond + 1; ++i) {
                    for (int j = 0; j < yDegreeSecond + 1; ++j) {
                        for (int h = 0; h < xDegreeFirst + 1; ++h) {
                            for (int k = 0; k < yDegreeFirst + 1; ++k) {
                                double val = 0.5 * referenceBoundaryRightLeft(k, j) * referenceIntegral(h, i);
                                localPenality(h * (yDegreeFirst + 1) + k, i * (yDegreeSecond + 1) + j) = val;
                            }
                        }
                    }

                }
                yDomainFirst = firstElement.getLengthY();
                yDomainSecond = secondElement.getLengthY();
                localPenality *= 4/sqrt(yDomainFirst)/sqrt(yDomainSecond)/yDomainSecond;
            break;

            default:
                throw std::invalid_argument("The accepted value are only 0, 1, 2 or 3.");
            break;
        }
    }
    return localPenality;
}
