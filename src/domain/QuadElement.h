//
// Created by Riccardo on 22/10/2017.
// Class for a 2D rectangular element of a mesh

#ifndef PACS_PELI_CHIAPPA_QUADELEMENT_H
#define PACS_PELI_CHIAPPA_QUADELEMENT_H


#include <tuple>
#include "Edge.h"
#include <vector>

class quadElement {
private:
    std::tuple<int,int,int,int> indices_;
    std::vector<edge> edges_;

public:
    quadElement();
    quadElement(std::tuple<int,int,int,int> indices,std::vector<edge> edges);
    quadElement(std::tuple<int,int,int,int> indices);
    ~quadElement();

    std::tuple<int,int,int,int> getIndices() const;
    void setIndices(std::tuple<int,int,int,int> const indices);

    std::vector<edge> getEdges() const;
    edge getEdge(int index) const ;
    void setEdges(std::vector<edge> const edges);
    void setEdgeNormal(const int index, Eigen::Vector2d normal );
    void setEdgeAdjacent(const int index, int const adjacentElem);

};

#endif //PACS_PELI_CHIAPPA_QUADELEMENT_H
