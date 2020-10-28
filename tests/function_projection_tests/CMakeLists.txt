 include_directories("${EIGEN3_INCLUDE_DIR}" lib/googletest-master/googletest/include
         lib/googletest-master/googlemock/include lib/googletest-master/googletest)

 add_executable(runFunctionProjectionTest ../../src/basis/HierarchicalBasis.cpp FunctionProjectionTest.cpp
         ../../src/utilities/FunctionProjection.h ../../src/utilities/FunctionProjection.cpp ../../src/utilities/GaussLobattoQuadrature.cpp
         ../../src/utilities/GaussLobattoQuadrature.h ../../src/basis/HierarchicalBasis.cpp ../../src/basis/HierarchicalBasis.h
         ../../src/basis/Basis.h ../../src/basis/Basis.cpp)

 target_link_libraries(runFunctionProjectionTest gtest gtest_main)
