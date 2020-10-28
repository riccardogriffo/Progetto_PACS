//
// Created by Alberto Chiappa on 19/09/17.
//

#include "Basis.h"


Basis::Basis() = default;

Eigen::MatrixXd Basis::evaluate(Eigen::VectorXd &csi) {
    return Eigen::MatrixXd::Zero(csi.size(), 1);
}

std::pair<double, double> Basis::getDomain() {
    return std::make_pair(-1., 1.);
}


