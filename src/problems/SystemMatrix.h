//
// Created by Riccardo on 04/01/2018.
//

#ifndef PACS_PELI_CHIAPPA_SYSTEMMATRIX_H
#define PACS_PELI_CHIAPPA_SYSTEMMATRIX_H

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include<vector>
#include "../domain/FeSpace.h"

typedef Eigen::Triplet<double> T;

class SystemMatrix {
private:
    Eigen::SparseMatrix<double> matrix_;
    std::vector<T> triplets_;
    Eigen::VectorXi blockStartingPoints_;
    Eigen::VectorXi blockDimensions_;
public:
    explicit SystemMatrix(FeSpace feSpace);
    void addBlock(const Eigen::MatrixXd &block,const int &indexFirstElem, const int &indexSecondElem);
    Eigen::SparseMatrix<double> getMatrix();
    void buildSparseMatrix();







};


#endif //PACS_PELI_CHIAPPA_SYSTEMMATRIX_H
