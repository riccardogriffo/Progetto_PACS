//
// Created by Alberto Chiappa on 15/12/17.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../../src/utilities/FunctionProjection.h"
#include "../../src/basis/HierarchicalBasis.h"

namespace {
    class FunctionProjectionTest: public testing::Test {
    public:
        std::function<double (double, double)> xy;

        std::shared_ptr<Basis> xBasis;

        std::shared_ptr<Basis> yBasis;

        FunctionProjectionTest(): xy([] (double x, double y) -> double { return x*y; }),
                                  xBasis(new HierarchicalBasis(3)),
                                  yBasis(new HierarchicalBasis(3)) {};
    };
}

TEST_F(FunctionProjectionTest, testProjectsConstantCorrectly) {
    auto constant = [] (double, double) -> double {return 1;};

    FunctionProjection functionProjection(constant, xBasis, yBasis);

    Eigen::VectorXd coeff = functionProjection.computeProjection(3);

    ASSERT_NEAR(coeff(1), 0., 1e-10);

}