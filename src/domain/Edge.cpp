//
// Created by Riccardo on 04/11/2017.
//

#include "Edge.h"

edge::edge() {

}

edge::edge(int startingNodeIndex, int endingNodeIndex, Eigen::Vector2d normal):
        startingNodeIndex_(startingNodeIndex),
        endingNodeIndex_(endingNodeIndex),
        normal_(normal){};


edge::edge(int startingNodeIndex, int endingNodeIndex):startingNodeIndex_(startingNodeIndex),endingNodeIndex_(endingNodeIndex){};

edge::~edge()=default;

int edge::getStartingNode() const{
    return startingNodeIndex_;
}

int edge::getEndingNode() const{
    return endingNodeIndex_;
}

Eigen::Vector2d edge::getNormal() const{
    return normal_;
}

int edge::getAdjacent() const {
    return adjacentElem_;
}

void edge::setStartingNode(const int &startingNodeIndex) {
    startingNodeIndex_=startingNodeIndex;
}

void edge::setEndingNode(const int &endingNodeIndex) {
    endingNodeIndex_=endingNodeIndex;
}

void edge::setNormal(const Eigen::Vector2d &normal) {
    normal_=normal;
}

void edge::setAdjacentElem(const int adjacentElem) {
    adjacentElem_=adjacentElem;
}

double edge::lenght(const node &startingNode, const node &endingNode) const {
    return startingNode.distance(endingNode);
}

Eigen::Vector2d edge::computeNormal(const node &startingNode, const node &endingNode) {
    Eigen::Vector2d normal;
    normal << endingNode.getY()-startingNode.getY(), startingNode.getX()-endingNode.getX();
    normal.normalize();
    return normal;

}

bool edge::hasAdjacentElem() {
    return adjacentElem_ == -1 ? false : true;
}










