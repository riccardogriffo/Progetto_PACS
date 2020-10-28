//
// Created by Alberto Chiappa on 25/09/17.
//

#include <iostream>
#include "GradientHierarchicalBasis.h"
#include "HierarchicalBasis.h"

#define EPS 1e-10


GradientHierarchicalBasis::GradientHierarchicalBasis(unsigned int numFunctions,
                                                     std::pair<double, double> domain):
        numFunctions(numFunctions), domain(domain) {}


Eigen::MatrixXd GradientHierarchicalBasis::evaluate(Eigen::VectorXd &csi) {
/*!   -------------------------------------------------------------------------
!   Computes the gradient of hierarchical (Legendre) basis functions on the points
!   csi.
!
!   Note: 1) if csi is of dimension M, then PHI and gradPHI  must be of dimension
!            MxNr, otherwise an error message is given and the program stopped.
!         2) PHI(:,k), i.e. the k-th column of PHI contains the
!            values of the k-th Legendre polinomial on the M points csi(0:M-1)
!         3) PHI(:,k), i.e. the k-th column of gradPHI contains the
!            values of the derivative of the k-th Legendre polynomial on the M
!            points csi(0:M-1)
!-------------------------------------------------------------------------*/

    if (csi.minCoeff() <= domain.first - EPS || csi.maxCoeff() >= domain.second + EPS) {
        throw std::invalid_argument("The points in which the functions are being evaluated exceed domain limits.");
    }

    Eigen::VectorXd csiShift = 2. / (domain.second - domain.first) *
                               (csi - (domain.second + domain.first) / 2. * Eigen::VectorXd::Ones(csi.size()));

    HierarchicalBasis hierarchicalBasis(numFunctions, std::make_pair(-1.,1.));

    Eigen::MatrixXd phi = hierarchicalBasis.evaluate(csiShift);

    Eigen::MatrixXd gradPhi(phi.rows(), phi.cols());

    gradPhi.col(0).setZero();
    if (numFunctions > 1) { // also other terms are needed
        gradPhi.col(1).setOnes();
        if (numFunctions > 2) { // also other terms are needed
            for (int i = 2; i < numFunctions; i++) {

                gradPhi.col(i) = 1./i * (sqrt(2*(2*i - 1)) * phi.col(i-1) + (2*i - 1) * csiShift.cwiseProduct(gradPhi.col(i-1)) -
                                  ((i-1) * gradPhi.col(i-2)));
            }
        }
    }

    for (int i = 0; i < numFunctions; i++) {
        gradPhi.col(i) *= sqrt((2*i+1.)/(domain.second-domain.first))*2/(domain.second-domain.first);
    }

    return gradPhi;
}

std::pair<double, double> GradientHierarchicalBasis::getDomain() {
    return this->domain;
}

