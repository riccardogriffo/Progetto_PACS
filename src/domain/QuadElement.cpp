//
// Created by Riccardo on 22/10/2017.
//

#include "QuadElement.h"
using namespace std;

quadElement::quadElement()=default;

quadElement::quadElement(tuple<int,int,int,int> indices, vector<edge> edges):indices_(indices),edges_(edges){};

quadElement::quadElement(std::tuple<int, int, int, int> indices):indices_(indices) {
};

quadElement::~quadElement()=default;

tuple<int,int,int,int> quadElement::getIndices() const{
    return indices_;
};

void quadElement::setIndices(tuple<int,int,int,int> const indices){
    indices_=indices;
}

std::vector<edge> quadElement::getEdges() const{
    return edges_;
};

void quadElement::setEdges(std::vector<edge> const edges){
    edges_=edges;
}

void quadElement::setEdgeNormal(const int index, Eigen::Vector2d normal) {
    edges_[index].setNormal(normal);
}

void quadElement::setEdgeAdjacent(const int index, const int adjacent) {
    edges_[index].setAdjacentElem(adjacent);
}

edge quadElement::getEdge(int index) const {
    return edges_[index];
}


