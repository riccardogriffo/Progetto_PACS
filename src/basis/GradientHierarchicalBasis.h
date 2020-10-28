//
// Created by Alberto Chiappa on 25/09/17.
//

#ifndef PACS_PELI_CHIAPPA_GRADIENTHIERARCHICALBASIS_H
#define PACS_PELI_CHIAPPA_GRADIENTHIERARCHICALBASIS_H


#include "Basis.h"

class GradientHierarchicalBasis: public Basis {
private:
    unsigned int numFunctions;

    const std::pair<double, double> domain;

public:
    explicit GradientHierarchicalBasis(unsigned int numFunctions,
                                       std::pair<double, double> domain = std::make_pair(-1., 1.));

    Eigen::MatrixXd evaluate(Eigen::VectorXd &csi) override;

    std::pair<double, double> getDomain() override;


};


#endif //PACS_PELI_CHIAPPA_GRADIENTHIERARCHICALBASIS_H
