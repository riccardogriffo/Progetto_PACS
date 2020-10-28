//
// Created by Riccardo on 01/12/2017.
//

#include "Operator.h"

Operator::Operator()=default;

void Operator::setReferenceElement(FeSpaceElement &referenceElement) {
    referenceElement_=referenceElement;
}

