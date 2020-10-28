//
// Created by Alberto Chiappa on 25/09/17.
//

#include <iostream>
#include "GradientLagrangeBasis.h"


#define EPS 1e-10


GradientLagrangeBasis::GradientLagrangeBasis(Eigen::VectorXd const &nodes): nodes(nodes) {}

Eigen::MatrixXd GradientLagrangeBasis::evaluate(Eigen::VectorXd &csi) {
    /*!   -------------------------------------------------------------------------
!   Computes the gradients (first derivative) of the  Lagrange basis functions.
!   Remark: if M is the number of points on which basis functions are to be
!   evaluated, and Nr is the number of basis functions, then, in principle, these
!   integers M and Nr could be different, and in fact this is the case for
!   the evaluation of the basis functions on the enriched grid
!   csi(M) = vector of points on which basis function are to be evaluated
!   csi_nod(Nr) = vector of interpolation nodes of the basis
!   ------------------------------------------------------------------------------*/
    // Initialize the matrix with the evaluation of the basis functions

    double firstNode = nodes.minCoeff();
    double lastNode = nodes.maxCoeff();

    if (csi.minCoeff() <= firstNode - EPS || csi.maxCoeff() >= lastNode + EPS) {
        throw std::invalid_argument("The points in which the functions are being evaluated exceed domain limits.");
    }

    Eigen::VectorXd csiShift = 2/(lastNode-firstNode) *
                               (csi - (lastNode+firstNode)/2*Eigen::VectorXd::Ones(csi.size()));
    Eigen::VectorXd nodesShift = 2/(lastNode-firstNode) *
                                 (nodes - (lastNode+firstNode)/2*Eigen::VectorXd::Ones(nodes.size()));

    Eigen::MatrixXd gradPhi = Eigen::MatrixXd::Ones(csi.size(), nodes.size());

    Eigen::VectorXd prod(csi.size());

    gradPhi.setZero();  //initialization to zero

    for(int i = 0; i < nodes.size(); i++) { // cycle over basis functions
        for(int s = 0; s < nodes.size(); s++) { // cycle over the nodes
            if(s != i) {
                prod.setOnes(); // initialization (1 is the neutral element for the multiplication)
                for(int j = 0; j < nodes.size(); j++) { // internal cycle over the nodes to compute the product of non-derived facors
                    if ((j != i) && (j != s)) {
                        prod = (prod.cwiseProduct(csiShift) - nodesShift(j) * prod) / (nodesShift(i) - nodesShift(j));
                    }
                }
                gradPhi.col(i) = gradPhi.col(i) + 1./(nodesShift(i)-nodesShift(s)) * prod * 2 / (lastNode - firstNode);

            }
        }
    }

    return gradPhi;
}


std::pair<double, double> GradientLagrangeBasis::getDomain() {
    return std::make_pair(nodes.minCoeff(), nodes.maxCoeff());
}

