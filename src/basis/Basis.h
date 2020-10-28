//
// Created by Alberto Chiappa on 19/09/17.
//

#ifndef PACS_PELI_CHIAPPA_BASIS_H
#define PACS_PELI_CHIAPPA_BASIS_H


#include <Eigen/Dense>

class Basis {

public:
    Basis();

    virtual Eigen::MatrixXd evaluate(Eigen::VectorXd &csi);

    virtual std::pair<double, double> getDomain();

};
#endif //PACS_PELI_CHIAPPA_BASIS_H
