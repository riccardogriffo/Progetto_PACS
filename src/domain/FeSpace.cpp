//
// Created by Riccardo on 12/11/2017.
//

#include "FeSpace.h"

FeSpace::FeSpace() {

}

FeSpace::FeSpace(unsigned int maxDegree, quadMesh mesh,
                 Eigen::Matrix<int, 2, Eigen::Dynamic> elemDegrees):
    maxDegree_(maxDegree),
    mesh_(mesh){
    feSpaceElements_=std::vector<FeSpaceElement>(mesh.getNumberElements());
    numberBasisFunctions_={0};
    for (int i=0; i<mesh.getNumberElements();++i){
        feSpaceElements_[i]=FeSpaceElement(elemDegrees(0,i),elemDegrees(1,i), mesh_.getXDomain(mesh.getElement(i)), mesh_.getYDomain(mesh.getElement(i)));
        numberBasisFunctions_+= (elemDegrees(0,i)+1)*(elemDegrees(1,i)+1);
    }
}

FeSpace::~FeSpace() =default;

quadMesh FeSpace::getMesh() {
    return mesh_;
}

unsigned int FeSpace::getMaxDegree() const{
    return maxDegree_;
}

std::vector<FeSpaceElement> FeSpace::getFeSpaceElements() {
    return feSpaceElements_;
}

int FeSpace::getNumberBasisFunctions() {
    return numberBasisFunctions_;
}





