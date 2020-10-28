 include_directories("${EIGEN3_INCLUDE_DIR}" lib/googletest-master/googletest/include
         lib/googletest-master/googlemock/include lib/googletest-master/googletest)

 add_executable(runSimpleTest SimpleTest.cpp)

 target_link_libraries(runSimpleTest gtest gtest_main)
