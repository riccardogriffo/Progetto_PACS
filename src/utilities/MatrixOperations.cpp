//
// Created by Alberto Chiappa on 08/11/17.
//

#include "MatrixOperations.h"

Eigen::VectorXd MatrixOperations::vectorElementWiseProduct(Eigen::VectorXd first, Eigen::VectorXd second) {
        long numElements = first.size();
        if(second.size() != numElements){
            throw std::invalid_argument("The two vectors must have the same size");
        }

        Eigen::VectorXd termProduct(numElements);

        for(int i = 0; i < numElements; i++) {
            termProduct(i) = first(i)*second(i);
        }
        return termProduct;
    }

Eigen::MatrixXd MatrixOperations::tensorProduct(const Eigen::MatrixXd &first, const Eigen::MatrixXd &second) {
    Eigen::MatrixXd result(first.rows()*second.rows(), first.cols()*second.cols());

    for (int i = 0; i < first.rows(); ++i){
        for (int j = 0; j < first.cols(); ++j) {
            for (int k = 0; k < second.rows(); ++k) {
                for (int h = 0; h < second.cols(); ++h) {
                    double value = first(i, j)*second(k, h);
                    long idx1 = first.rows()*k + i;
                    long idx2 = first.cols()*h + j;
                    result(idx1, idx2) = value;
                }
            }
        }
    }
    return result;
}
