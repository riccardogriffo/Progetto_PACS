//
// Created by Alberto Chiappa on 26/11/17.
//

#include <gtest/gtest.h>
#include <Eigen/Dense>
#include "../../src/basis/LagrangeBasis.h"
#include "../../src/basis/GradientLagrangeBasis.h"
#include "../../src/basis/HierarchicalBasis.h"
#include "../../src/basis/GradientHierarchicalBasis.h"

#define NUM_BASIS_FUNCTIONS 6
#define NUM_NODES 4
#define NUM_EVAL_POINTS 5
#define TOL 1e-10
#define LEFT 0.
#define RIGHT 1.

namespace {
    class SystemMatrixTest: public testing::Test {
    public:
        LagrangeBasis lagrangeBasis;
        GradientLagrangeBasis gradientLagrangeBasis;
        HierarchicalBasis hierarchicalBasis;
        GradientHierarchicalBasis gradientHierarchicalBasis;
        std::pair<double, double> domain;

        SystemMatrixTest(): domain(LEFT, RIGHT),
                     lagrangeBasis(Eigen::VectorXd::LinSpaced(NUM_NODES, LEFT, RIGHT)),
                     gradientLagrangeBasis(Eigen::VectorXd::LinSpaced(NUM_NODES, LEFT, RIGHT)),
                     hierarchicalBasis(NUM_BASIS_FUNCTIONS, std::make_pair(LEFT, RIGHT)),
                     gradientHierarchicalBasis(NUM_BASIS_FUNCTIONS, std::make_pair(LEFT, RIGHT)) {}
    };
}