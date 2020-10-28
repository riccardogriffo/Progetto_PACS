//
// Created by Riccardo on 22/10/2017.
//

#include "Node.h"
#include <cmath>
using namespace std;

node::node(){}

node::node(tuple<double,double> coordinates):coordinates_(coordinates){};

node::~node(){}

tuple<double,double> node::getCoord() const{
    return coordinates_;
}

double node::getX() const{
    return get<0>(coordinates_);
}

double node::getY() const{
    return get<1>(coordinates_);
}

void node::setCoord(tuple<double,double> coordinates){
    coordinates_=coordinates;
}

void node::setX(double x){
    get<0>(coordinates_)=x;
}

void node::setY(double y){
    get<0>(coordinates_ ) = y;
}

double node::distance(const node node2) const {
    return sqrt(get<0>(this->coordinates_)*node2.getX()+get<1>(this->coordinates_)*node2.getY());
}
