#include <gtest/gtest.h>
#include "storage/array2d.h"

class ArrayTest : public testing::Test
{
protected:
    ArrayTest()
    {
        std::array<int, 2> size = {3, 3};
        Array2D arr(size);
    }

    Array2D arr;
};

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