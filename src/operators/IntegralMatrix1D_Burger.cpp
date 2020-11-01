//
// Created by Alberto Chiappa on 24/11/17.
//

#include "IntegralMatrix1D.h"
#include "../utilities/MatrixOperations.h"

IntegralMatrix1D::IntegralMatrix1D(std::shared_ptr<Basis> firstBasisPtr, std::shared_ptr<Basis> secondBasisPtr):
        firstBasisPtr(std::move(firstBasisPtr)), secondBasisPtr(std::move(secondBasisPtr)) {

    if (this->firstBasisPtr->getDomain() != this->secondBasisPtr->getDomain()) {
        throw std::invalid_argument("The two bases have to be defined on the same domain.");
    }

    this->quad = GaussLobattoQuadrature();
}

Eigen::MatrixXd IntegralMatrix1D::generate(unsigned int maxPolynomialDegree) {
    Eigen::VectorXd quadNodes;
    Eigen::VectorXd quadWeights;

    unsigned int numNodes = (maxPolynomialDegree + 3)/2 + ((maxPolynomialDegree + 3) % 2 != 0);
    std::tie(quadNodes, quadWeights) = quad.computeNodesAndWeights(numNodes, firstBasisPtr->getDomain());

    Eigen::MatrixXd polyMatrix1 = firstBasisPtr->evaluate(quadNodes);
    Eigen::MatrixXd polyMatrix2 = secondBasisPtr->evaluate(quadNodes);
    Eigen::MatrixXd result = Eigen::MatrixXd::Ones(polyMatrix1.cols(), polyMatrix2.cols());
    for (int i = 0; i < result.rows(); i++) {
        for (int j = 0; j < result.cols(); j++) {
            result(i, j) = MatrixOperations().vectorElementWiseProduct(polyMatrix1.col(i), polyMatrix2.col(j)).dot(quadWeights);
        }
    }
    return result;
}




