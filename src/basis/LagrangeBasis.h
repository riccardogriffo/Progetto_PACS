//
// Created by Alberto Chiappa on 20/09/17.
//

#ifndef PACS_PELI_CHIAPPA_LAGRANGEBASIS_H
#define PACS_PELI_CHIAPPA_LAGRANGEBASIS_H


#include "Basis.h"

class LagrangeBasis: public Basis {
private:
    Eigen::VectorXd const nodes;

public:

    explicit LagrangeBasis(Eigen::VectorXd const &nodes);

    ~LagrangeBasis();

    Eigen::MatrixXd evaluate(Eigen::VectorXd &csi) override;

    std::pair<double, double> getDomain() override;

};


#endif //PACS_PELI_CHIAPPA_LAGRANGEBASIS_H
