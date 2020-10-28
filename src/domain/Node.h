//
// Created by Riccardo on 22/10/2017.
// Class of a 2D node of a mesh

#ifndef PACS_PELI_CHIAPPA_NODE_H
#define PACS_PELI_CHIAPPA_NODE_H

#include <tuple>


class node {
private:
    std::tuple<double, double> coordinates_;
public:
    node();
    node(std::tuple<double,double> coordinates);
    ~node();

    std::tuple<double,double> getCoord()const;
    double getX()const;
    double getY()const;
    void setCoord(std::tuple<double,double> coordinates);
    void setX(double x);
    void setY(double y);
    double distance(const node node2) const;



};

#endif //PACS_PELI_CHIAPPA_NODE_H
