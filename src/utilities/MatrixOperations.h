//
// Created by Alberto Chiappa on 08/11/17.
//

#ifndef PACS_PELI_CHIAPPA_MATRIXVECTOROPERATIONS_H
#define PACS_PELI_CHIAPPA_MATRIXVECTOROPERATIONS_H


#include <Eigen/Dense>
#include <stdexcept>

class MatrixOperations {
public:
    static Eigen::VectorXd vectorElementWiseProduct(Eigen::VectorXd first, Eigen::VectorXd second);

    static Eigen::MatrixXd tensorProduct(const Eigen::MatrixXd &first, const Eigen::MatrixXd &second);

};


#endif //PACS_PELI_CHIAPPA_MATRIXVECTOROPERATIONS_H
