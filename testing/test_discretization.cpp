#include <gtest/gtest.h>
#include <array>
#include "storage/centralDifferences.h"

class StaggeredGridTest : public testing::Test
{
public:
    StaggeredGridTest() : grid({6, 5}, {1.0, 1.0})
    {
    }

protected:
    StaggeredGrid grid;
};

class CentralDifferencesTest : public testing::Test
{
public:
    CentralDifferencesTest() : diff({6, 5}, {1.0, 1.0})
    {
        // fill u
        int iBegin = diff.uIBegin();
        int iEnd = diff.uIEnd();
        int jBegin = diff.uJBegin();
        int jEnd = diff.uJEnd();
        for (int i = iBegin - 1)
        {
            for
        }

        // fill v

        // fill p
    }

protected:
    CentralDifferences diff;
};

TEST_F(StaggeredGridTest, TestUI)
{
    int iBegin = grid.uIBegin();
    int iEnd = grid.uIEnd();

    EXPECT_EQ(iBegin, 1);
    EXPECT_EQ(iEnd, 4);
}

TEST_F(StaggeredGridTest, TestVI)
{
    int iBegin = grid.vIBegin();
    int iEnd = grid.vIEnd();

    EXPECT_EQ(iBegin, 1);
    EXPECT_EQ(iEnd, 5);
}

TEST_F(StaggeredGridTest, TestUJ)
{
    int jBegin = grid.uJBegin();
    int jEnd = grid.uJEnd();

    EXPECT_EQ(jBegin, 1);
    EXPECT_EQ(jEnd, 4);
}

TEST_F(StaggeredGridTest, TestVJ)
{
    int jBegin = grid.vJBegin();
    int jEnd = grid.vJEnd();

    EXPECT_EQ(jBegin, 1);
    EXPECT_EQ(jEnd, 3);
}

TEST_F(StaggeredGridTest, TestPI)
{
    int iBegin = grid.pIBegin();
    int iEnd = grid.pIEnd();

    EXPECT_EQ(iBegin, 1);
    EXPECT_EQ(iEnd, 5);
}

TEST_F(StaggeredGridTest, TestPJ)
{
    int jBegin = grid.pJBegin();
    int jEnd = grid.pJEnd();

    EXPECT_EQ(jBegin, 1);
    EXPECT_EQ(jEnd, 4);
}

// class CentralDifferencesTest : public testing::Test
// {
// public:
//     CentralDifferencesTest() : disc(std::array<int, 2>{4, 4}, std::array<double, 2>{1.0, 1.0})
//     {
//         }

// protected:
//     CentralDifferences disc;
// };

// TEST_F(CentralDifferencesTest, TestDu2Dx)
// {
// }