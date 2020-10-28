//
// Created by Alberto Chiappa on 23/10/17.
//

#include <gtest/gtest.h>
#include "ClassName.h"

namespace {
    class ClassDeclaration: public testing::Test {
    public:
        ClassName obj;

        ClassDeclaration() = default;
    };
}

TEST_F(ClassDeclaration, nameOfTheFirstTest) {
    ASSERT_EQ("", "");
}

TEST_F(ClassDeclaration, nameOfSecondTest) {
    int a = 1;
    ASSERT_EQ(1, a);
}

TEST_F(ClassDeclaration, testField1Initialization){
    ASSERT_EQ(obj.getField1(), 1);
}
