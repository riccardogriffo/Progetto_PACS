//
// Created by Alberto Chiappa on 25/09/17.
//

#ifndef PACS_PELI_CHIAPPA_GRADIENTLAGRANGEBASIS_H
#define PACS_PELI_CHIAPPA_GRADIENTLAGRANGEBASIS_H


#include "Basis.h"

class GradientLagrangeBasis: public Basis {
private:
    Eigen::VectorXd const nodes;
public:

    explicit GradientLagrangeBasis(Eigen::VectorXd const &nodes);

    Eigen::MatrixXd evaluate(Eigen::VectorXd &csi) override;

    std::pair<double, double> getDomain() override;

};


#endif //PACS_PELI_CHIAPPA_GRADIENTLAGRANGEBASIS_H
