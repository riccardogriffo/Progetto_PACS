//
// Created by Riccardo on 01/12/2017.
//

#ifndef PACS_PELI_CHIAPPA_OPERATOR_H
#define PACS_PELI_CHIAPPA_OPERATOR_H

#include"../basis/Basis.h"
#include"../domain/FeSpaceElement.h"
#include <Eigen/Dense>
#include <utility>
#include <memory>
#include<iostream>


class Operator {
protected:
    Eigen::MatrixXd referenceMatrix1D_;
    FeSpaceElement referenceElement_;
    //FeSpace fespace_;

public:
    Operator();
    void setReferenceElement(FeSpaceElement &referenceElement);
    virtual void computeReferenceMatrix1D() = 0;
    // virtual void updateSystemMatrix(SystemMatrix &systemMatrix)=0;
    //virtual Eigen::MatrixXd computeLocalMatrix(FeSpaceElement feSpaceElement);
};


#endif //PACS_PELI_CHIAPPA_OPERATOR_H
