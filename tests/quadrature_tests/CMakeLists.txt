 include_directories("${EIGEN3_INCLUDE_DIR}" lib/googletest-master/googletest/include
         lib/googletest-master/googlemock/include lib/googletest-master/googletest)

 add_executable(runQuadratureTest QuadratureTest.cpp ../../src/utilities/GaussLobattoQuadrature.cpp
         ../../src/utilities/GaussLobattoQuadrature.h ../../src/utilities/MatrixOperations.h ../../src/utilities/MatrixOperations.cpp)

 target_link_libraries(runQuadratureTest gtest gtest_main)
