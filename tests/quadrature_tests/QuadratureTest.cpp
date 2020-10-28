//
// Created by Alberto Chiappa on 26/11/17.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../../src/utilities/GaussLobattoQuadrature.h"

#define TOL 1e-10
#define LEFT 0.
#define RIGHT 1.
#define ACCURACY_DEGREE 10

namespace {
    class QuadratureTest: public testing::Test {
    public:
        std::pair<double, double> domain;
        GaussLobattoQuadrature quad;
        Eigen::VectorXd quadNodes;
        Eigen::VectorXd quadWeights;
        unsigned int numNodes;

        QuadratureTest(): domain(LEFT, RIGHT), quad() {

            numNodes = (ACCURACY_DEGREE + 3)/2 + ((ACCURACY_DEGREE + 3) % 2 != 0);

            std::tie(quadNodes, quadWeights) = quad.computeNodesAndWeights(numNodes, domain);
        }
    };

    TEST_F(QuadratureTest, testWeightsHaveRightSize) {
        ASSERT_EQ(quadWeights.size(), numNodes);
    }

    TEST_F(QuadratureTest, testNodesHaveRightSize) {
        ASSERT_EQ(quadNodes.size(), numNodes);
    }

    TEST_F(QuadratureTest, testIntegratesLinearFunctionCorrectly) {
        std::function<double (double)> linFunc = [](double x) { return x; };
        std::function<double (double)> quadFunc = [](double x) { return x*x; };

        double exactResult = 0.5*(quadFunc(RIGHT)-quadFunc(LEFT));
        double numericalResult = quad.integrate1d(linFunc, domain);

        ASSERT_NEAR(exactResult, numericalResult, TOL);
    }

    TEST_F(QuadratureTest, testIntegratesQuadraticFunctionCorrectly) {
        std::function<double (double)> quadFunc = [](double x) { return x*x; };
        std::function<double (double)> cubFunc = [](double x) { return x*x*x; };
        double exactResult = 1./3*(cubFunc(RIGHT)-cubFunc(LEFT));
        double numericalResult = quad.integrate1d(quadFunc, domain);

        ASSERT_NEAR(exactResult, numericalResult, TOL);
    }

    TEST_F(QuadratureTest, testIntegrateConstant2D) {
        std::function<double (double, double)> constFunc = [](double x, double y) {return 1.;};

        double exactResult = (RIGHT-LEFT)*(RIGHT-LEFT);
        double numericalResult = quad.integrate2d(constFunc, domain, domain);

        ASSERT_NEAR(exactResult, numericalResult, TOL);
    }

    TEST_F(QuadratureTest, testIntegrateXY2D) {
        std::function<double (double, double)> xy = [](double x, double y) {return x*y;};
        std::function<double (double)> quadFunc = [](double x) { return x*x; };

        double exactResult = 1./4*(quadFunc(RIGHT)-quadFunc(LEFT))*(quadFunc(RIGHT)-quadFunc(LEFT));
        double numericalResult = quad.integrate2d(xy, domain, domain);

        ASSERT_NEAR(exactResult, numericalResult, TOL);
    }
}