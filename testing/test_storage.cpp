#include <gtest/gtest.h>
#include <array>
#include "storage/array2d.h"
#include "storage/fieldVariable.h"

class ArrayTest : public testing::Test
{
public:
    ArrayTest() :
        arr({3, 3})
    {
        // std::array<int,2> size{3, 3};
        // Array2D arr = Array2D::Array2D(size);
    }

protected:
    Array2D arr;
};

class FieldVariableTest : public testing::Test
{
protected:
    FieldVariableTest()
    {
        std::array<int, 2> size = {3,3};
        std::array<double, 2> originU = {0, 0};
        std::array<double, 2> originP = {0, -0.25};
        std::array<double, 2> meshWidthU = {1, 1};
        std::array<double, 2> meshWidthP = {0.5, 0.5};

        varU = std::make_unique<FieldVariable>(size, originU, meshWidthU);
        varP = std::make_unique<FieldVariable>(size, originP, meshWidthP);

        for (int i = 0; i < size[0]; i++)
        {
            for (int j = 0; j < size[1]; j++)
            {
                (*varU)(i, j) = i;
                (*varP)(i, j) = std::max(i, j);
            }
        }
    }

    std::unique_ptr<FieldVariable> varU;
    std::unique_ptr<FieldVariable> varP;
    
}; 

TEST(Arraytest, CreateArray)
{
    std::array<int,2> size{3,3};
    Array2D arr = Array2D(size);

}

TEST_F(ArrayTest, GetSize)
{
    std::array<int, 2> expected = {3, 3};
    std::array<int, 2> size = arr.size();

    EXPECT_EQ(expected, size);
}

TEST_F(ArrayTest, SetGetValue)
{
    arr(2, 2) = 50;
    arr(0, 2) = 20;

    double val1 = arr(2, 2);
    double val2 = arr(0, 2);

    EXPECT_EQ(val1, 50);
    EXPECT_EQ(val2, 20);
}

TEST_F(ArrayTest, SetValueOurofBounds)
{
    EXPECT_DEATH({ arr(5, 2) = 200; }, "");
}

TEST_F(ArrayTest, GetValueOutOfBounds)
{
    EXPECT_DEATH({ double val = arr(5, 2); }, "");
}

TEST_F(FieldVariableTest, interpolateUAtCorner)
{
    double val = varU->interpolateAt(0, 3);
    EXPECT_EQ(val, 0);
    
}

TEST_F(FieldVariableTest, interpolatePAtCorner)
{
    double val = varP->interpolateAt(0, 3 * 0.5 - 0.25);
    EXPECT_EQ(val, 0);
}


TEST_F(FieldVariableTest, interpolateUAtMiddle)
{
    double val = varU->interpolateAt(1.5, 1.5);
    EXPECT_EQ(val, 1);
}

TEST_F(FieldVariableTest, interpolatePAtMiddle)
{
    double val = varP->interpolateAt(0.75, 0.75);
    EXPECT_EQ(val, 1.5);
}

TEST_F(FieldVariableTest, interpolateUOutOfBounds)
{
    EXPECT_DEATH({varP->interpolateAt(-0.5, -0.5);}, "");

}

TEST_F(FieldVariableTest, interpolatePOutOfBounds)
{
    EXPECT_DEATH({varU->interpolateAt(-0.5, -0.5);}, "");
}
