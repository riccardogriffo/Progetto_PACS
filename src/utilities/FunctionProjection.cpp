//
// Created by Alberto Chiappa on 05/12/17.
//

#include "FunctionProjection.h"


FunctionProjection::FunctionProjection(const std::function<double(double, double)> &func2d,
                                       std::shared_ptr<Basis> xBasisPtr,
                                       std::shared_ptr<Basis> yBasisPtr): func2d(func2d),
                                                                          xBasisPtr(std::move(xBasisPtr)),
                                                                          yBasisPtr(std::move(yBasisPtr)){
    quad = GaussLobattoQuadrature();
}

Eigen::VectorXd FunctionProjection::computeProjection(unsigned int accuracyDegree) {

    unsigned int numNodes = (accuracyDegree + 3)/2 + ((accuracyDegree) % 2 == 0);


    Eigen::VectorXd xQuadNodes;
    Eigen::VectorXd yQuadNodes;
    Eigen::VectorXd xQuadWeights;
    Eigen::VectorXd yQuadWeights;

    std::tie(xQuadNodes, xQuadWeights) = quad.computeNodesAndWeights(numNodes, xBasisPtr->getDomain());
    std::tie(yQuadNodes, yQuadWeights) = quad.computeNodesAndWeights(numNodes, yBasisPtr->getDomain());

    Eigen::MatrixXd xPolyBasis = xBasisPtr->evaluate(xQuadNodes);
    Eigen::MatrixXd yPolyBasis = yBasisPtr->evaluate(yQuadNodes);

    Eigen::VectorXd projCoeff=Eigen::VectorXd::Zero(xPolyBasis.cols()*yPolyBasis.cols());

    for (int i = 0; i < xPolyBasis.cols(); ++i) {
        for (int j = 0; j < yPolyBasis.cols(); ++j) {
            for (int k = 0; k < xPolyBasis.rows(); ++k) {
                for (int h = 0; h < yPolyBasis.rows(); ++h) {

                    double partialsum = func2d(xQuadNodes(k), yQuadNodes(h))*xPolyBasis(k, i)*yPolyBasis(h, j)*xQuadWeights(k)*yQuadWeights(h);
                    projCoeff(i*yPolyBasis.cols() + j) += partialsum;
                }
            }
        }
    }

    return projCoeff;
}

