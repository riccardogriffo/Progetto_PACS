//
// Created by Riccardo on 12/11/2017.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include<tuple>
#include "../../src/domain/QuadMesh.h"


namespace {
    class MeshTest : public testing::Test {
    public:
        quadMesh mesh;
        int nX = 3;
        int nY = 2;
        double hX = 1;
        double hY = 1;

        MeshTest() {

            mesh = quadMesh(nX, nY, hX, hY);
        };
    };
}


TEST_F(MeshTest, checkNumberElements) {
    ASSERT_EQ(mesh.getElements().size(), (nX-1)*(nY-1));
}

TEST_F(MeshTest, checkNumberNodes) {
    ASSERT_EQ(mesh.getNodes().size(), nX*nY);
}


TEST_F(MeshTest, checkNormals) {
    double eps=1e-4;
    auto elements=mesh.getElements();
    for (auto it=elements.begin(); it<elements.end();++it){
        auto edges=it->getEdges();
        ASSERT_NEAR(edges[0].getNormal()[0],0,eps);
        ASSERT_NEAR(edges[0].getNormal()[1],-1,eps);
        ASSERT_NEAR(edges[1].getNormal()[0],1,eps);
        ASSERT_NEAR(edges[1].getNormal()[1],0,eps);
        ASSERT_NEAR(edges[2].getNormal()[0],0,eps);
        ASSERT_NEAR(edges[2].getNormal()[1],1,eps);
        ASSERT_NEAR(edges[3].getNormal()[0],-1,eps);
        ASSERT_NEAR(edges[3].getNormal()[1],0,eps);
    }
}

TEST_F(MeshTest, checkEdges) {
    auto elements=mesh.getElements();
    auto edges=elements[0].getEdges();
    ASSERT_EQ(edges[0].getStartingNode(),0);
    ASSERT_EQ(edges[0].getEndingNode(),2);
    ASSERT_EQ(edges[1].getStartingNode(),2);
    ASSERT_EQ(edges[1].getEndingNode(),3);
    ASSERT_EQ(edges[2].getStartingNode(),3);
    ASSERT_EQ(edges[2].getEndingNode(),1);
    ASSERT_EQ(edges[3].getStartingNode(),1);
    ASSERT_EQ(edges[3].getEndingNode(),0);

    edges=elements[1].getEdges();
    ASSERT_EQ(edges[0].getStartingNode(),2);
    ASSERT_EQ(edges[0].getEndingNode(),4);
    ASSERT_EQ(edges[1].getStartingNode(),4);
    ASSERT_EQ(edges[1].getEndingNode(),5);
    ASSERT_EQ(edges[2].getStartingNode(),5);
    ASSERT_EQ(edges[2].getEndingNode(),3);
    ASSERT_EQ(edges[3].getStartingNode(),3);
    ASSERT_EQ(edges[3].getEndingNode(),2);

}

TEST_F(MeshTest, checkAdjacent) {
    auto elements=mesh.getElements();
    auto edges=elements[0].getEdges();
    ASSERT_EQ(edges[0].getAdjacent(),-1);
    ASSERT_EQ(edges[1].getAdjacent(),1);
    ASSERT_EQ(edges[2].getAdjacent(),-1);
    ASSERT_EQ(edges[3].getAdjacent(),-1);

    edges=elements[1].getEdges();
    ASSERT_EQ(edges[0].getAdjacent(),-1);
    ASSERT_EQ(edges[1].getAdjacent(),-1);
    ASSERT_EQ(edges[2].getAdjacent(),-1);
    ASSERT_EQ(edges[3].getAdjacent(),0);

}

TEST_F(MeshTest, ZeroROws){
    nX=0;
    nY=4;
    hX=1;
    hY=1;
    quadMesh zeromesh(nX,nY,hX,hY);
    ASSERT_EQ(zeromesh.getElements().size(),0);
    ASSERT_EQ(zeromesh.getNodes().size(),nY);


}