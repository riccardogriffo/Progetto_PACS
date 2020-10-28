//
// Created by Alberto Chiappa on 20/09/17.
//

#include <iostream>
#include "LagrangeBasis.h"

#define EPS 1e-10

LagrangeBasis::LagrangeBasis(Eigen::VectorXd const &nodes): nodes(nodes) {}

Eigen::MatrixXd LagrangeBasis::evaluate(Eigen::VectorXd &csi) {
    double firstNode = nodes.minCoeff();
    double lastNode = nodes.maxCoeff();

    if (csi.minCoeff() <= firstNode - EPS || csi.maxCoeff() >= lastNode + EPS) {
        throw std::invalid_argument("The points in which the functions are being evaluated exceed domain limits.");
    }

    Eigen::VectorXd csiShift = 2/(lastNode-firstNode) *
                               (csi - (lastNode+firstNode)/2*Eigen::VectorXd::Ones(csi.size()));
    Eigen::VectorXd nodesShift = 2/(lastNode-firstNode) *
                               (nodes - (lastNode+firstNode)/2*Eigen::VectorXd::Ones(nodes.size()));

    // Initialize the matrix with the evaluation of the basis functions
    Eigen::MatrixXd phi = Eigen::MatrixXd::Ones(csi.size(), nodes.size());

    // Basis functions are computed according to Lagrange formula

    for(int i = 0; i < nodes.size(); i++) {   // cycle over basis functions
        for (int j = 0; j < nodes.size(); j++) { // cycle over the nodes
            if(j != i)
                phi.col(i) = (phi.col(i).cwiseProduct(csiShift)-nodesShift(j)*phi.col(i))/(nodesShift(i)-nodesShift(j));


        }
    }

    return phi;
}

std::pair<double, double> LagrangeBasis::getDomain() {
    return std::make_pair(nodes.minCoeff(), nodes.maxCoeff());
}

LagrangeBasis::~LagrangeBasis() = default;
