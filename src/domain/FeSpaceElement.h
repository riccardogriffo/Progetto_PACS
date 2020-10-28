//
// Created by Riccardo on 01/12/2017.
//

#ifndef PACS_PELI_CHIAPPA_REFERENCEELEMENT_H
#define PACS_PELI_CHIAPPA_REFERENCEELEMENT_H

#include <Eigen/Dense>
#include "QuadElement.h"
#include "../basis/HierarchicalBasis.h"
#include "../basis/GradientHierarchicalBasis.h"
#include <math.h>

class FeSpaceElement {
private:
    int degreeX_;
    int degreeY_;
    std::pair<double ,double> domainX_;
    std::pair<double ,double> domainY_;

public:
    FeSpaceElement();
    FeSpaceElement(int degreeX,int degreeY,std::pair<double ,double> domainX,std::pair<double ,double> domainY );
    ~FeSpaceElement();

    int getDegreeX() const;
    int getDegreeY() const;
    std::pair<double ,double> getDomainX() const;
    std::pair<double ,double> getDomainY() const;

    double getLengthX()const;
    double getLengthY()const;

    Eigen::MatrixXd getBasisFunctionsX(Eigen::VectorXd evalPoints);
    Eigen::MatrixXd getBasisFunctionsY(Eigen::VectorXd evalPoints);

    Eigen::MatrixXd getGradBasisFunctionsX(Eigen::VectorXd evalPoints);
    Eigen::MatrixXd getGradBasisFunctionsY(Eigen::VectorXd evalPoints);

};


#endif //PACS_PELI_CHIAPPA_REFERENCEELEMENT_H
