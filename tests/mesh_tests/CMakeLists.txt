 include_directories("${EIGEN3_INCLUDE_DIR}" lib/googletest-master/googletest/include
         lib/googletest-master/googlemock/include lib/googletest-master/googletest)

 add_executable(runMeshTest MeshTest.cpp ../../src/domain/QuadMesh.cpp ../../src/domain/Node.cpp ../../src/domain/Edge.cpp ../../src/domain/QuadElement.cpp)

 target_link_libraries(runMeshTest gtest gtest_main)
