//
// Created by Alberto Chiappa on 12/12/17.
//

#ifndef PACS_PELI_CHIAPPA_RIGHTHANDSIDEFUNCTION_H
#define PACS_PELI_CHIAPPA_RIGHTHANDSIDEFUNCTION_H


#include <functional>
#include "../domain/FeSpace.h"
#include <Eigen/Dense>


class RightHandSideFunction {
private:
    std::function<double (double, double)> func2d;
    FeSpace fespace;
    unsigned int accuracyDegree;

public:
    RightHandSideFunction(FeSpace fespace,
                          std::function<double (double, double)> func2d,
                          unsigned int accuracyDegree = 10);
    Eigen::VectorXd getSystemVector();

};


#endif //PACS_PELI_CHIAPPA_RIGHTHANDSIDEFUNCTION_H
