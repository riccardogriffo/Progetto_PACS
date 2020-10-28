//
// Created by Alberto Chiappa on 25/09/17.
//

#include <iostream>
#include <stdexcept>
#include <math.h>
#include "HierarchicalBasis.h"

#define EPS 1e-10

HierarchicalBasis::HierarchicalBasis(unsigned int numFunctions,
                                     std::pair<double, double> domain):
        numFunctions(numFunctions), domain(domain){}

Eigen::MatrixXd HierarchicalBasis::evaluate(Eigen::VectorXd &csi) {
    //  -------------------------------------------------------------------------
    //   Computes the hierarchical (Legendre) basis functions on the points csi.
    //
    //   Note: 1) if csi is of dimension M, then PHI must be of dimension MxNr, other-
    //           wise an error message is given and the program stopped.
    //         2) PHI(:,k), i.e. the k-th column of PHI contains the
    //           values of the k-th Legendre polinomial on the M points csi(0:M-1)
    //   -------------------------------------------------------------------------

    if (csi.minCoeff() <= domain.first - EPS || csi.maxCoeff() >= domain.second + EPS) {
        throw std::invalid_argument("The points in which the functions are being evaluated exceed domain limits.");
    }

    Eigen::VectorXd csiShift = 2/(domain.second-domain.first) *
            (csi - (domain.second+domain.first)/2*Eigen::VectorXd::Ones(csi.size()));
    Eigen::MatrixXd phi(csi.size(), numFunctions);


    phi.col(0).setOnes(); // primo polinomio di Legendre

    if (numFunctions > 1) {    //also other terms are needed
        phi.col(1) = csiShift;
        // Secondo polinomio di Legendre
        if(numFunctions > 2) {
            // also other terms are needed
            for (int i = 2; i < numFunctions; i++) { // three terms recurrence formula
                phi.col(i) = 1./i * ((2*i-1) * csiShift.cwiseProduct(phi.col(i-1)) - (i-1) * phi.col(i-2));
            }
        }
    }

    for (int i = 0; i < numFunctions; i++) {
        phi.col(i) *= sqrt((2*i+1.)/(domain.second-domain.first));
    }

    return phi;
}

std::pair<double, double> HierarchicalBasis::getDomain() {
    return this->domain;
}

