//
// Created by Riccardo on 22/10/2017.
//Class of a 2D mesh, composed of rectangles

#include "QuadMesh.h"
#include "Edge.h"
#include <iostream>
using namespace std;

quadMesh::quadMesh() =default;

//constructor for uniform mesh of rectangle
quadMesh::quadMesh(int nX, int nY, double hX, double hY) {
    if(nX>0 && nY>0)
        nodes_= vector<node>(nX*nY);
    else if(nX>0)
        nodes_=vector<node>(nX);
    else if(nY>0)
        nodes_=vector<node>(nY);

    if(nX>0 && nY>0){
        elements_= vector<quadElement>((nX-1)*(nY-1));
    }
    else
        elements_=vector<quadElement>(0);
    nX_=nX;
    nY_=nY;
    hX_=hX;
    hY_=hY;

    domainX_ = std::make_pair(0., hX*(nX-1));
    domainY_ = std::make_pair(0., hY*(nY-1));

    for (int i=0;i<nX;++i){
        for(int j=0;j<nY;++j){
            nodes_[i*nY+j]=node(make_tuple(hX*i,hY*j));
        }
    }
    for (int i=0;i<nX-1;++i) {
        for (int j = 0; j < nY - 1; ++j) {
            elements_[i * (nY - 1) + j] = quadElement(
                    make_tuple(i * nY + j, (i + 1) * nY + j, (i + 1) * nY + j + 1, i * nY + j + 1));
        }
    }
    computeEdges();
    computeNormals();
    computeAdjacents();
}

quadMesh::quadMesh(std::pair<double,double> domainX,
                   std::pair<double ,double> domainY,
                   int nX,
                   int nY): domainX_(domainX), domainY_(domainY){
    double hX= ((domainX.second-domainX.first)/(nX-1));
    double hY= ((domainY.second-domainY.first)/(nY-1));
    if(nX>0 && nY>0)
        nodes_= vector<node>(nX*nY);
    else if(nX>0)
        nodes_=vector<node>(nX);
    else if(nY>0)
        nodes_=vector<node>(nY);

    if(nX>0 && nY>0){
        elements_= vector<quadElement>((nX-1)*(nY-1));
    }
    else
        elements_=vector<quadElement>(0);
    nX_=nX;
    nY_=nY;
    hX_=hX;
    hY_=hY;
    for (int i=0;i<nX;++i){
        for(int j=0;j<nY;++j){
            nodes_[i*nY+j]=node(make_tuple(domainX.first+hX*i,domainY.first+ hY*j));
        }
    }
    for (int i=0;i<nX-1;++i) {
        for (int j = 0; j < nY - 1; ++j) {
            elements_[i * (nY - 1) + j] = quadElement(
                    make_tuple(i * nY + j, (i + 1) * nY + j, (i + 1) * nY + j + 1, i * nY + j + 1));
        }
    }
    computeEdges();
    computeNormals();
    computeAdjacents();
}

quadMesh::~quadMesh() {}

std::vector<node> quadMesh::getNodes() const {
    return nodes_;
}

node quadMesh::getNode(const int & index) {
    return nodes_[index];
}

std::vector<quadElement> quadMesh::getElements() const {
    return elements_;
}

quadElement quadMesh::getElement(const int & index) {
    return elements_[index];
}

int quadMesh::getNX() const {
    return nX_;
}

int quadMesh::getNY() const {
    return nY_;
}

void quadMesh::setNodes(const std::vector<node> & nodes) {
    nodes_=nodes;
}

void quadMesh::setElements(const std::vector<quadElement> & elements) {
    elements_=elements;
}

void quadMesh::setNX(const int & nX) {
    nX_=nX;
}

void quadMesh::setNY(const int & nY) {
    nY_=nY;
}

bool quadMesh::isBorder(const int index) {
    for (int i=0;i<nX_;++i){
        if(index==i*nY_ || index==i*nY_+nY_-1)
            return true;
    }
    for(int j=1; j<nY_-1;++j){
        if(index==j || index==(nX_-1)*nY_+j)
            return true;
    }
    return false;
}

bool quadMesh::isBorder(const edge edge, const int i) {
    if(nX_ >2 && nY_ >2){
       return  (isBorder(edge.getStartingNode()) && isBorder(edge.getEndingNode()));
    }
    else if(nY_==2){
        if (i==0 || i==2)
            return true;
        else
            return ((edge.getStartingNode()==0 && edge.getEndingNode()==1) || (edge.getStartingNode()==nX_*nY_-2 && edge.getEndingNode()==nX_*nY_-1));
    }
    else if(nX_==2){
        if(i==1 || i==3)
            return true;
        else
            return ((edge.getStartingNode()==0 && edge.getEndingNode()==nY_) || (edge.getStartingNode()==nX_*nY_-1 && edge.getEndingNode()==nY_-1));
    }
    else{
        return true;
    }
}

double quadMesh::computeLength(const edge &edge) {
    return edge.lenght(nodes_[edge.getStartingNode()], nodes_[edge.getEndingNode()]);
}

unsigned int quadMesh::getNumberNodes() {
    return nodes_.size();
}

unsigned int quadMesh::getNumberElements() {
    return elements_.size();
}

std::pair<double, double> quadMesh::getXDomain(const quadElement &elem) {
    auto indices=elem.getIndices();
    return pair<double, double>(nodes_[get<0>(indices)].getX(),nodes_[get<1>(indices)].getX());
}

std::pair<double, double> quadMesh::getYDomain(const quadElement &elem) {
    auto indices=elem.getIndices();
    return pair<double, double>(nodes_[get<1>(indices)].getY(),nodes_[get<2>(indices)].getY());
}


double quadMesh::elementArea( const quadElement &elem) {
    auto edges=elem.getEdges();
    return computeLength(edges[0])+computeLength(edges[1]);
}

void quadMesh::printInfo(){
    cout<<"Rectangular mesh"<<endl<<"Number of vertices: "<<nodes_.size()<<endl<<"Number of rectangles: "<<elements_.size()<<endl;
    cout<<"Nodes:"<<endl;
    printNodes();
    cout<<"Elements:"<<endl;
    printElements();
}

void quadMesh::printNodes() {
    for (int i=0; i<nodes_.size();++i){
        cout<<i<<" "<<nodes_[i].getX()<<" "<<nodes_[i].getY()<<endl;
    }
}

void quadMesh::printElements() {
    for (int i=0; i<elements_.size();++i){
        tuple<int,int,int,int> indices=elements_[i].getIndices();
        cout<<get<0>(indices)<<" "<<get<1>(indices)<<" "<<get<2>(indices)<<" "<<get<3>(indices)<<endl;
        auto edges=elements_[i].getEdges();
        for(int j=0;j<edges.size();++j){
            cout<<"edge "<<edges[j].getStartingNode()<<" "<<edges[j].getEndingNode()<<endl;
            cout<<"normal "<<endl<<edges[j].getNormal()<<endl;
            cout <<"adjacent "<<edges[j].getAdjacent()<<endl;
        }
    }
}

void quadMesh::printBorder() {
    for (int i=0; i<nodes_.size();++i){
        if(isBorder(i))
            cout<<i<<endl;
    }
}
void quadMesh::computeEdges() {
    for (auto it=elements_.begin();it<elements_.end(); ++it){
        auto indices=it->getIndices();
        std::vector<edge> edges={edge(get<0>(indices),get<1>(indices)), edge(get<1>(indices),get<2>(indices)),
                                 edge(get<2>(indices),get<3>(indices)), edge(get<3>(indices),get<0>(indices))};
        it->setEdges(edges);
    }
}

void quadMesh::computeNormals() {
    for (auto it=elements_.begin();it<elements_.end(); ++it){
        auto edges = it->getEdges();
        for (int i=0; i<edges.size(); i++){
            auto normal=edges[i].computeNormal(nodes_[edges[i].getStartingNode()],nodes_[edges[i].getEndingNode()]);
            it->setEdgeNormal(i,normal);
        }
    }
}

void quadMesh::computeAdjacents() {
    for (int j=0; j<elements_.size();++j){
        auto edges=elements_[j].getEdges();
        elements_[j].setEdgeAdjacent(0,j-1);
        elements_[j].setEdgeAdjacent(1,j+nY_-1);
        elements_[j].setEdgeAdjacent(2,j+1);
        elements_[j].setEdgeAdjacent(3,j-nY_+1);
        for (int i=0; i<edges.size();++i) {
            if (isBorder(edges[i],i)) {
                elements_[j].setEdgeAdjacent(i, -1);
            }
        }
    }

}

std::pair<double, double> quadMesh::getDomainX() {
    return domainX_;
}

std::pair<double, double> quadMesh::getDomainY() {
    return domainY_;
}
