//
// Created by Riccardo on 22/10/2017.
// Class of a 2D mesh, composed of rectangles

#ifndef PACS_PELI_CHIAPPA_QUADMESH_H
#define PACS_PELI_CHIAPPA_QUADMESH_H

#include <vector>
#include "Node.h"
#include "QuadElement.h"
#include "FeSpaceElement.h"


class quadMesh {
private:
    std::vector<node> nodes_;
    std::vector<quadElement> elements_;
    int nX_; //number of vertices along the x axis
    int nY_; //number of vertices along the y axis
    double hX_;
    double hY_;
    std::pair<double,double> domainX_;
    std::pair<double ,double> domainY_;

public:
    quadMesh();
    quadMesh(int nX, int nY, double hX, double hY);
    quadMesh(std::pair<double,double> domainX, std::pair<double ,double> domainY, int nX, int nY);
    ~quadMesh();

    std::vector<node> getNodes() const;
    node getNode(const int &index);
    std::vector<quadElement> getElements() const;
    quadElement getElement(const int &index);
    int getNX() const;
    int getNY() const;
    std::pair<double,double> getDomainX();
    std::pair<double,double> getDomainY();

    void setNodes(const std::vector<node> & nodes);
    void setElements(const std::vector<quadElement> & elements);
    void setNX(const int & nRows);
    void setNY(const int & nCols);

    std::pair<double,double> getXDomain(const quadElement &elem);
    std::pair<double,double> getYDomain(const quadElement &elem);

    bool isBorder(const int index);
    bool isBorder(const edge edge, const int index); //arguments: edge and index of the edge in the element(0,1,2,3)
    double computeLength(const edge &edge);
    double elementArea(const quadElement &elem);

    unsigned int getNumberElements();
    unsigned int getNumberNodes();


    void computeEdges();
    void computeNormals();
    void computeAdjacents();

    void printInfo();
    void printNodes();
    void printElements();
    void printBorder();
};


#endif //PACS_PELI_CHIAPPA_QUADMESH_H
