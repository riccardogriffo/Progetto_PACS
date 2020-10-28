//
// Created by Riccardo on 04/01/2018.
//

#include "SystemMatrix.h"
#define EPS 1e-10
SystemMatrix::SystemMatrix(FeSpace feSpace) {

    auto feSpaceElemVec = feSpace.getFeSpaceElements();
    blockDimensions_ = Eigen::VectorXi::Zero(feSpaceElemVec.size());
    blockStartingPoints_ = Eigen::VectorXi::Zero(feSpaceElemVec.size());
    int startingIndex = 0;
    for (int i = 0; i < feSpaceElemVec.size(); ++i) {
        blockStartingPoints_[i] = startingIndex;
        blockDimensions_[i] = (feSpaceElemVec[i].getDegreeX() + 1) * (feSpaceElemVec[i].getDegreeY() + 1);
        startingIndex += blockDimensions_[i];
    }
    matrix_ = Eigen::SparseMatrix<double,Eigen::RowMajor>(feSpace.getNumberBasisFunctions(), feSpace.getNumberBasisFunctions());
    int estimation_of_entries = feSpaceElemVec.size()*5*(feSpace.getMaxDegree()+1)*(feSpace.getMaxDegree()+1); //CORREGGERE!!!!
    triplets_.reserve(estimation_of_entries);


}

void SystemMatrix::addBlock(const Eigen::MatrixXd &block, const int &indexFirstElem, const int &indexSecondElem) {

    for (int i =0; i<block.rows();++i){
        for(int j=0;j<block.cols();++j){
            if(abs(block(i,j))>EPS) {
                triplets_.push_back(
                        T(blockStartingPoints_[indexFirstElem] + i, blockStartingPoints_[indexSecondElem] + j,
                          block(i, j)));
            }
        }
    }
}

Eigen::SparseMatrix<double> SystemMatrix::getMatrix() {
    return matrix_;
}

void SystemMatrix::buildSparseMatrix(){
    matrix_.setFromTriplets(triplets_.begin(), triplets_.end());
}

