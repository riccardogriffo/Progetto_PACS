//
// Created by Riccardo on 04/11/2017.
//

#ifndef PACS_PELI_CHIAPPA_EDGE_H
#define PACS_PELI_CHIAPPA_EDGE_H

#include <Eigen/Dense>
#include "Node.h"

class edge {
private:
    int startingNodeIndex_;
    int endingNodeIndex_;
    Eigen::Vector2d normal_; //outward-pointing normal (with respect to the quadElement)
    int adjacentElem_;  //index of the adjacent rectangle
public:
    edge();
    edge(int startingNodeIndex, int endingNodeIndex, Eigen::Vector2d normal);
    edge(int startingNodeIndex, int endingNodeIndex);
    ~edge();

    int getStartingNode()const;
    int getEndingNode()const;
    Eigen::Vector2d getNormal()const;
    int getAdjacent() const;
    void setStartingNode(const int &startingNodeIndex);
    void setEndingNode(const int &endingNodeIndex);
    void setNormal(const Eigen::Vector2d &normal);
    void setAdjacentElem(const int adjacentElem);
    double lenght(const node &startingNode, const node &endingNode) const;
    Eigen::Vector2d computeNormal(const node &startingNode, const node &endingNode);
    bool hasAdjacentElem();


};


#endif //PACS_PELI_CHIAPPA_EDGE_H
