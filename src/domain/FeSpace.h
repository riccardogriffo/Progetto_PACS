//
// Created by Riccardo on 12/11/2017.
//

#ifndef PACS_PELI_CHIAPPA_FESPACE_H
#define PACS_PELI_CHIAPPA_FESPACE_H

#include"QuadMesh.h"
#include "FeSpaceElement.h"
#include <string>
#include <Eigen/Dense>
#include <iostream>

class FeSpace {
private:
    unsigned int maxDegree_; //global maximum degree of the basis
    quadMesh mesh_;
    std::vector<FeSpaceElement> feSpaceElements_;
    int numberBasisFunctions_;
public:
    FeSpace();
    FeSpace(unsigned int maxDegree, quadMesh mesh, Eigen::Matrix <int, 2, Eigen::Dynamic> elemDegrees_);
    ~FeSpace();
    unsigned int getMaxDegree() const;
    quadMesh getMesh();
    std::vector<FeSpaceElement> getFeSpaceElements();
    int getNumberBasisFunctions();

};


#endif //PACS_PELI_CHIAPPA_FESPACE_H
