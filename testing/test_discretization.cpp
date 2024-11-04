#include <gtest/gtest.h>
#include <array>
#include "storage/centralDifferences.h"
#include <math.h>

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

        for (int i = iBegin - 1; i <= iEnd; i++)
        {
            for (int j = jBegin - 1; j <= jEnd; j++)
            {
                diff.u(i, j) = std::max(i, j);
            }
        }

        // fill v
        iBegin = diff.vIBegin();
        iEnd = diff.vIEnd();
        jBegin = diff.vJBegin();
        jEnd = diff.vJEnd();

        for (int i = iBegin - 1; i <= iEnd; i++)
        {
            for (int j = jBegin - 1; j <= jEnd; j++)
            {
                diff.v(i, j) = std::max(i, j);
            }
        }

        // fill p
        iBegin = diff.pIBegin();
        iEnd = diff.pIEnd();
        jBegin = diff.pJBegin();
        jEnd = diff.pJEnd();

        for (int i = iBegin - 1; i <= iEnd; i++)
        {
            for (int j = jBegin - 1; j <= jEnd; j++)
            {
                diff.p(i, j) = std::max(i, j);
            }
        }
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

TEST_F(CentralDifferencesTest, Du2DxTest)
{
    int iBegin = diff.uIBegin();
    int iEnd = diff.uIEnd();

    // left border
    double left = diff.computeDu2Dx(iBegin, 0);
    // right border
    double right = diff.computeDu2Dx(iEnd, 0);
    // middle
    double middle = diff.computeDu2Dx(2, 0);

    EXPECT_EQ(left, 2.0);
    EXPECT_EQ(right, std::pow(2, 2) - std::pow(3.5, 2));
    EXPECT_EQ(middle, std::pow(2.5, 2) - std::pow(1.5, 2));
}

TEST_F(CentralDifferencesTest, Dv2DyTest)
{
    int jBegin = diff.vJBegin();
    int jEnd = diff.vJEnd();

    // bottom border
    double bottom = diff.computeDv2Dy(2, jBegin);
    // top border
    double top = diff.computeDv2Dy(2, jEnd);

    EXPECT_EQ(bottom, 0);
    EXPECT_EQ(top, std::pow(1.5, 2) - std::pow(2.5, 2));
}

TEST_F(CentralDifferencesTest, DuvDxTest)
{
}

TEST_F(CentralDifferencesTest, DuvDyTest)
{
}

TEST_F(CentralDifferencesTest, D2uDx2Test)
{
}

TEST_F(CentralDifferencesTest, D2uDy2Test)
{
}

TEST_F(CentralDifferencesTest, D2vDx2Test)
{
}

TEST_F(CentralDifferencesTest, D2vDy2Test)
{
}

TEST_F(CentralDifferencesTest, DpDxTest)
{
}

TEST_F(CentralDifferencesTest, DpDyTest)
{
}
