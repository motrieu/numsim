#include "donorCell.h"
#include <cmath>

DonorCell::DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha) :
    Discretization(nCells, meshWidth), alpha_(alpha)
{
}

double DonorCell::computeDu2Dx(int i, int j) const
{
    const double uRightMean = (u(i,j) + u(i+1,j)) / 2.0;
    const double uLeftMean = (u(i-1,j) + u(i,j)) / 2.0;
    const double uRightDifference = (u(i,j) - u(i+1,j)) / 2.0;
    const double uLeftDifference = (u(i-1,j) - u(i,j)) / 2.0;

    return (uRightMean*uRightMean - uLeftMean*uLeftMean) / dx()
                + (alpha_/dx()) * (std::fabs(uRightMean)*uRightDifference - std::fabs(uLeftMean)*uLeftDifference);
   
}
 	
double DonorCell::computeDv2Dy(int i, int j) const
{
    const double vUpperMean = (v(i,j) + v(i,j+1)) / 2.0;
    const double vLowerMean = (v(i,j-1) + v(i,j)) / 2.0;
    const double vUpperDifference = (v(i,j) - v(i,j+1)) / 2.0;
    const double vLowerDifference = (v(i,j-1) - v(i,j)) / 2.0;

    return (vUpperMean*vUpperMean - vLowerMean*vLowerMean) / dy()
                + (alpha_/dy()) * (std::fabs(vUpperMean)*vUpperDifference - std::fabs(vLowerMean)*vLowerDifference);
}

double DonorCell::computeDuvDx(int i, int j) const
{
    const double uRightMean = (u(i,j) + u(i,j+1)) / 2.0;
    const double vRightMean = (v(i,j) + v(i+1,j)) / 2.0;
    const double uLeftMean = (u(i-1,j) + u(i-1,j+1)) / 2.0;
    const double vLeftMean = (v(i-1,j) + v(i,j)) / 2.0;

    const double vRightDifference = (v(i,j) - v(i+1,j)) / 2.0;
    const double vLeftDifference = (v(i-1,j) - v(i,j)) / 2.0;

    return (uRightMean*vRightMean - uLeftMean*vLeftMean) / dx()
                + (alpha_/dx()) * (std::fabs(uRightMean)*vRightDifference - std::fabs(uLeftMean)*vLeftDifference);
}

double DonorCell::computeDuvDy(int i, int j) const
{
    const double vUpperMean = (v(i,j) + v(i+1,j)) / 2.0;
    const double uUpperMean = (u(i,j) + u(i,j+1)) / 2.0;
    const double vLowerMean = (v(i,j-1) + v(i+1,j-1)) / 2.0;
    const double uLowerMean = (u(i,j-1) + u(i,j)) / 2.0;

    const double uUpperDifference = (u(i,j) - u(i,j+1)) / 2.0;
    const double uLowerDifference = (u(i,j-1) - u(i,j)) / 2.0;
    
    return (vUpperMean*uUpperMean - vLowerMean*uLowerMean) / dy()
            + (alpha_/dy()) * (std::fabs(vUpperMean)*uUpperDifference - std::fabs(vLowerMean)*uLowerDifference);
}