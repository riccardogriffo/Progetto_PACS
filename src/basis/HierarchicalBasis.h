//
// Created by Alberto Chiappa on 25/09/17.
//

#ifndef PACS_PELI_CHIAPPA_HIERARCHICALBASIS_H
#define PACS_PELI_CHIAPPA_HIERARCHICALBASIS_H


#include "Basis.h"

class HierarchicalBasis: public Basis {
private:
    const unsigned int numFunctions;
    const std::pair<double, double> domain;

public:
    explicit HierarchicalBasis(unsigned int numFunctions, std::pair<double, double> domain = std::make_pair(-1., 1.));

    Eigen::MatrixXd evaluate(Eigen::VectorXd &csi) override;

    std::pair<double, double> getDomain() override;

};


#endif //PACS_PELI_CHIAPPA_HIERARCHICALBASIS_H
