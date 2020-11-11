cmake_minimum_required(VERSION 3.8)
project(pacs)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-O2")


 ## SET( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
 SET( EIGEN3_INCLUDE_DIR "/usr/local/include/eigen3" )


 IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable ##EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
 ENDIF()


if (APPLE)
   add_definitions(-DGTEST_USE_OWN_TR1_TUPLE)
    add_definitions(-D__GLIBCXX__)
endif (APPLE)

 ##  add_subdirectory(tests)
 include_directories("${EIGEN3_INCLUDE_DIR}"
##        lib/googletest-master/googletest/include
##        lib/googletest-master/googlemock/include
          )

set(SOURCE_FILES main.cpp
                 src/basis/Basis.cpp
                 src/basis/Basis.h
                 src/basis/LagrangeBasis.cpp
                 src/basis/LagrangeBasis.h
                 src/basis/HierarchicalBasis.cpp
                 src/basis/HierarchicalBasis.h
                 src/basis/GradientLagrangeBasis.cpp
                 src/basis/GradientLagrangeBasis.h
                 src/basis/GradientHierarchicalBasis.cpp
                 src/basis/GradientHierarchicalBasis.h
                 src/utilities/GaussLobattoQuadrature.cpp
                 src/utilities/GaussLobattoQuadrature.h
                 src/utilities/MatrixOperations.cpp
                 src/utilities/MatrixOperations.h
                 src/domain/QuadElement.h
                 src/domain/QuadElement.cpp
                 src/domain/QuadMesh.h
                 src/domain/QuadMesh.cpp
                 src/domain/Node.h
                 src/domain/Node.cpp
                 src/domain/Edge.cpp
                 src/domain/Edge.h
                 src/domain/FeSpace.cpp
                 src/domain/FeSpace.h
                 src/operators/IntegralMatrix1D.cpp
                 src/operators/IntegralMatrix1D.h
                 src/utilities/FunctionProjection.cpp
                 src/utilities/FunctionProjection.h
 ##       tests/basis_tests/BasisTest.cpp tests/example_tests/ClassName.h tests/example_tests/SimpleTest.cpp
 ##       tests/system_matrix_tests/SystemMatrixTest.cpp
 ##       tests/quadrature_tests/QuadratureTest.cpp
 ##       tests/function_projection_tests/FunctionProjectionTest.cpp
                 src/operators/Operator.cpp
                 src/operators/Operator.h
                 src/domain/FeSpaceElement.cpp
                 src/domain/FeSpaceElement.h
                 src/operators/DiffusionOperator.cpp
                 src/operators/DiffusionOperator.h
                 src/operators/AdvectionOperator.cpp
                 src/operators/AdvectionOperator.h
                 src/operators/ReactionOperator.cpp
                 src/operators/ReactionOperator.h
                 src/operators/RightHandSideFunction.cpp
                 src/operators/RightHandSideFunction.h
##        src/operators/InteriorPenalityOperatorOld.cpp src/operators/InteriorPenalityOperatorOld.h
                 src/operators/RightHandSideFunction.cpp
                 src/operators/RightHandSideFunction.h
                 src/operators/StabilizerOperator.cpp
                 src/operators/StabilizerOperator.h
                 src/problems/SystemMatrix.cpp
                 src/problems/SystemMatrix.h
                 src/utilities/SolutionEvaluation.cpp
                 src/utilities/SolutionEvaluation.h
                 src/utilities/Error.cpp
                 src/utilities/Error.h
                 src/problems/SystemMatrix.cpp
                 src/problems/SystemMatrix.h
                 src/utilities/SolutionEvaluation.cpp
                 src/utilities/SolutionEvaluation.h
                 src/utilities/Error.cpp src/utilities/Error.h
                 src/operators/InteriorPenalityOperator.cpp
                 src/operators/InteriorPenalityOperator.h
                 src/operators/InteriorPenalitySymmetric.cpp
                 src/operators/InteriorPenalitySymmetric.h
                 src/problems/LaplaceProblem.cpp
                 src/problems/LaplaceProblem.h
                 src/problems/AdvectionDiffusionProblem.cpp
                 src/problems/AdvectionDiffusionProblem.h
                 src/operators/TransportEdgeOperator.h
                 src/operators/TransportEdgeOperator.cpp
                 )

add_executable(pacs ${SOURCE_FILES})

## target_link_libraries(pacs gtest gtest_main)
target_link_libraries(pacs )
