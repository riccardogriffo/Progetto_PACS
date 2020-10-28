 include_directories("${EIGEN3_INCLUDE_DIR}" lib/googletest-master/googletest/include
         lib/googletest-master/googlemock/include lib/googletest-master/googletest)

 add_executable(runBasisTest BasisTest.cpp ../../src/basis/LagrangeBasis.cpp ../../src/basis/LagrangeBasis.h
         ../../src/basis/Basis.h ../../src/basis/Basis.cpp ../../src/basis/HierarchicalBasis.cpp ../../src/basis/HierarchicalBasis.h
         ../../src/basis/GradientHierarchicalBasis.cpp ../../src/basis/GradientHierarchicalBasis.h
         ../../src/basis/GradientLagrangeBasis.cpp ../../src/basis/GradientLagrangeBasis.h ../../src/utilities/GaussLobattoQuadrature.cpp
         ../../src/utilities/GaussLobattoQuadrature.h ../../src/utilities/MatrixOperations.h ../../src/utilities/MatrixOperations.cpp)

 target_link_libraries(runBasisTest gtest gtest_main)
