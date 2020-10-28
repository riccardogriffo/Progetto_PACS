//
// Created by Alberto Chiappa on 11/11/17.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../../src/basis/LagrangeBasis.h"
#include "../../src/basis/GradientLagrangeBasis.h"
#include "../../src/basis/HierarchicalBasis.h"
#include "../../src/basis/GradientHierarchicalBasis.h"
#include "../../src/utilities/GaussLobattoQuadrature.h"
#include "../../src/utilities/MatrixOperations.h"


#define NUM_BASIS_FUNCTIONS 6
#define NUM_NODES 4
#define NUM_EVAL_POINTS 5
#define TOL 1e-10
#define LEFT -1.
#define RIGHT 2.

namespace {
    class BasisTest: public testing::Test {
    public:
        LagrangeBasis lagrangeBasis;
        GradientLagrangeBasis gradientLagrangeBasis;
        HierarchicalBasis hierarchicalBasis;
        GradientHierarchicalBasis gradientHierarchicalBasis;
        std::pair<double, double> domain;

        BasisTest(): domain(LEFT, RIGHT),
                     lagrangeBasis(Eigen::VectorXd::LinSpaced(NUM_NODES, LEFT, RIGHT)),
                     gradientLagrangeBasis(Eigen::VectorXd::LinSpaced(NUM_NODES, LEFT, RIGHT)),
                     hierarchicalBasis(NUM_BASIS_FUNCTIONS, std::make_pair(LEFT, RIGHT)),
                     gradientHierarchicalBasis(NUM_BASIS_FUNCTIONS, std::make_pair(LEFT, RIGHT)) {}
    };
}

TEST_F(BasisTest, testEvaluationLagrangeMatrixSize) {
    Eigen::VectorXd csi = Eigen::VectorXd::LinSpaced(NUM_EVAL_POINTS, domain.first, domain.second);
    Eigen::MatrixXd eval = lagrangeBasis.evaluate(csi);
    ASSERT_EQ(eval.rows(), NUM_EVAL_POINTS);
    ASSERT_EQ(eval.cols(), NUM_NODES);
}

TEST_F(BasisTest, testEvaluationGradLagrangeMatrixSize) {
    Eigen::VectorXd csi = Eigen::VectorXd::LinSpaced(NUM_EVAL_POINTS, domain.first, domain.second);
    Eigen::MatrixXd eval = gradientLagrangeBasis.evaluate(csi);
    ASSERT_EQ(eval.rows(), NUM_EVAL_POINTS);
    ASSERT_EQ(eval.cols(), NUM_NODES);
}

TEST_F(BasisTest, testEvaluationHierarchicalMatrixSize) {
    Eigen::VectorXd csi = Eigen::VectorXd::LinSpaced(NUM_EVAL_POINTS, domain.first, domain.second);
    Eigen::MatrixXd eval = hierarchicalBasis.evaluate(csi);
    std::cout << std::endl << eval << std::endl;
    ASSERT_EQ(eval.rows(), NUM_EVAL_POINTS);
    ASSERT_EQ(eval.cols(), NUM_BASIS_FUNCTIONS);
}

TEST_F(BasisTest, testEvaluationGradHierarchicalMatrixSize) {
    Eigen::VectorXd csi = Eigen::VectorXd::LinSpaced(NUM_EVAL_POINTS, domain.first, domain.second);
    Eigen::MatrixXd eval = gradientHierarchicalBasis.evaluate(csi);
    std::cout << std::endl << std::endl << eval << std::endl;
    ASSERT_EQ(eval.rows(), NUM_EVAL_POINTS);
    ASSERT_EQ(eval.cols(), NUM_BASIS_FUNCTIONS);
}

TEST_F(BasisTest, testLagrangeBasisIsZeroOrOneAtNodes) {
    Eigen::VectorXd csi = Eigen::VectorXd::LinSpaced(NUM_NODES, domain.first, domain.second);
    Eigen::MatrixXd eval = lagrangeBasis.evaluate(csi);

    for (int i = 0; i < NUM_NODES; ++i) {
        for (int j = 0; j < NUM_NODES; ++j) {
            if (i == j) {
                ASSERT_NEAR(eval(i, j), 1, TOL);
            }
            else {
                ASSERT_NEAR(eval(i, j), 0, TOL);
            }
        }
    }
}

TEST_F(BasisTest, testHierarchicalPolynomialsAreOrthogonal) {
    GaussLobattoQuadrature quad;
    Eigen::VectorXd quadNodes;
    Eigen::VectorXd quadWeights;

    std::tie(quadNodes, quadWeights) = quad.computeNodesAndWeights(NUM_BASIS_FUNCTIONS , domain);

    auto eval = hierarchicalBasis.evaluate(quadNodes);

    for (int i = 0; i < NUM_BASIS_FUNCTIONS; ++i) {
        for(int j = i + 1; j < NUM_BASIS_FUNCTIONS; ++j) {
            double result = MatrixOperations().vectorElementWiseProduct(eval.col(i), eval.col(j)).dot(quadWeights);
            ASSERT_NEAR(result, 0, TOL);
        }
    }
}

TEST_F(BasisTest, testHierarchicalPolynomialsAreNormalized) {
    GaussLobattoQuadrature quad;
    Eigen::VectorXd quadNodes;
    Eigen::VectorXd quadWeights;

    std::tie(quadNodes, quadWeights) = quad.computeNodesAndWeights(NUM_BASIS_FUNCTIONS + 2 , domain);

    auto eval = hierarchicalBasis.evaluate(quadNodes);

    for (int i = 0; i < NUM_BASIS_FUNCTIONS; ++i) {
        double result = MatrixOperations().vectorElementWiseProduct(eval.col(i), eval.col(i)).dot(quadWeights);
        ASSERT_NEAR(result, 1, TOL);
    }
}