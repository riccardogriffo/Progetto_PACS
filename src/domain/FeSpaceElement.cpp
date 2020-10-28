//
// Created by Riccardo on 01/12/2017.
//

#include "FeSpaceElement.h"

FeSpaceElement::FeSpaceElement()=default;

FeSpaceElement::FeSpaceElement(int degreeX, int degreeY, std::pair<double, double> domainX,
                                   std::pair<double, double> domainY):degreeX_(degreeX),degreeY_(degreeY),domainX_(domainX),domainY_(domainY) {}

FeSpaceElement::~FeSpaceElement() = default;

int FeSpaceElement::getDegreeX() const{
    return degreeX_;
}

int FeSpaceElement::getDegreeY() const {
    return degreeY_;
}

std::pair<double, double> FeSpaceElement::getDomainX() const {
    return domainX_;
}

std::pair<double, double> FeSpaceElement::getDomainY() const {
    return domainY_;
}

double FeSpaceElement::getLengthX() const{
    return std::abs(domainX_.first-domainX_.second);
}

double FeSpaceElement::getLengthY() const{
    return std::abs(domainY_.first-domainY_.second);
}

Eigen::MatrixXd FeSpaceElement::getBasisFunctionsX(Eigen::VectorXd evalPoints) {
    HierarchicalBasis basis(degreeX_ + 1, domainX_);
    return basis.evaluate(evalPoints);
}

Eigen::MatrixXd FeSpaceElement::getBasisFunctionsY(Eigen::VectorXd evalPoints) {
    HierarchicalBasis basis(degreeY_ + 1, domainY_);
    return basis.evaluate(evalPoints);
}

Eigen::MatrixXd FeSpaceElement::getGradBasisFunctionsX(Eigen::VectorXd evalPoints) {
    GradientHierarchicalBasis basis(degreeX_ + 1, domainX_);
    return basis.evaluate(evalPoints);
}

Eigen::MatrixXd FeSpaceElement::getGradBasisFunctionsY(Eigen::VectorXd evalPoints) {
    GradientHierarchicalBasis basis(degreeY_ + 1, domainY_);
    return basis.evaluate(evalPoints);
}





