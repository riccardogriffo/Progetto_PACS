//
// Created by Alberto Chiappa on 26/10/17.
//

#ifndef PACS_PELI_CHIAPPA_GAUSSLOBATTOQUADRATURE_H
#define PACS_PELI_CHIAPPA_GAUSSLOBATTOQUADRATURE_H


#include <Eigen/Dense>

class GaussLobattoQuadrature {
private:
    void LegendreCoefficients(unsigned int M, Eigen::VectorXd & a,Eigen::VectorXd & b);

public:
    explicit GaussLobattoQuadrature() = default;

    std::pair<Eigen::VectorXd, Eigen::VectorXd>  computeNodesAndWeights(unsigned int numNodes,
                                                                        std::pair<double, double> domain);

    template <typename Function>
    double integrate1d(Function integrand, std::pair<double, double> domain, unsigned int accuracyDegree = 10) {
        Eigen::VectorXd quadNodes;
        Eigen::VectorXd quadWeights;
        unsigned int numNodes = (accuracyDegree + 3)/2 + ((accuracyDegree + 3) % 2 != 0);

        std::tie(quadNodes, quadWeights) = computeNodesAndWeights(numNodes, domain);

        Eigen::VectorXd funcVals(numNodes);

        for (int i = 0; i < numNodes; ++i) {
            funcVals(i) = integrand(quadNodes(i));
        }

        return funcVals.dot(quadWeights);
    }

    template <typename Function2d>
    double integrate2d(Function2d integrand,
                       std::pair<double, double> xDomain,
                       std::pair<double, double> yDomain,
                       unsigned int xAccuracyDegree = 10,
                       unsigned int yAccuracyDegree = 10) {

        unsigned int xNumNodes = (xAccuracyDegree + 3)/2 + ((xAccuracyDegree + 3) % 2 != 0);
        unsigned int yNumNodes = (yAccuracyDegree + 3)/2 + ((yAccuracyDegree + 3) % 2 != 0);

        Eigen::VectorXd xQuadNodes;
        Eigen::VectorXd yQuadNodes;
        Eigen::VectorXd xQuadWeights;
        Eigen::VectorXd yQuadWeights;

        std::tie(xQuadNodes, xQuadWeights) = computeNodesAndWeights(xNumNodes, xDomain);
        std::tie(yQuadNodes, yQuadWeights) = computeNodesAndWeights(yNumNodes, yDomain);

        double result = 0;

        for (int i = 0; i < xNumNodes; ++i) {
            for (int j = 0; j < yNumNodes; ++j) {
                result += integrand(xQuadNodes(i), yQuadNodes(j)) * xQuadWeights(i) * yQuadWeights(j);
            }
        }

        return result;
    }
};


#endif //PACS_PELI_CHIAPPA_GAUSSLOBATTOQUADRATURE_H
