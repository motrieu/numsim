#include "donorCell.h"
#include <cmath>

double DonorCell::computeDu2Dx(int i, int j) const
{
    const double uRightMean = (u(i,j) + u(i+1,j)) / 2.0;
    const double uLeftMean = (u(i-1,j) + u(i,j)) / 2.0;
    const double uRightDifference = (u(i,j) - u(i+1,j)) / 2.0;
    const double uLeftDifference = (u(i-1,j) - u(i,j)) / 2.0;

    return (1/dx()) * (uRightMean*uRightMean - uLeftMean*uLeftMean)
                + (alpha_/dx()) * (std::abs(uRightMean)*uRightDifference - std::abs(uLeftMean)*uLeftDifference);
   
}
 	
double DonorCell::computeDv2Dy(int i, int j) const
{
    const double vUpperMean = (v(i,j) + v(i,j+1)) / 2.0;
    const double vLowerMean = (v(i,j-1) + v(i,j)) / 2.0;
    const double vUpperDifference = (v(i,j) - v(i,j+1)) / 2.0;
    const double vLowerDifference = (v(i,j-1) - v(i,j)) / 2.0;

    return (1.0/dy()) * (vUpperMean*vUpperMean - vLowerMean*vLowerMean) 
                + (alpha_/dy()) * (std::abs(vUpperMean)*vUpperDifference - std::abs(vLowerMean)*vLowerDifference);
}

double DonorCell::computeDuvDx(int i, int j) const
{
    const double uRightMean = (u(i,j) + u(i,j+1)) / 2.0;
    const double vRightMean = (v(i,j) + v(i+1,j)) / 2.0;
    const double uLeftMean = (u(i-1,j) + u(i-1,j+1)) / 2.0;
    const double vLeftMean = (v(i-1,j) + v(i,j)) / 2.0;

    const double vRightDifference = (v(i,j) - v(i+1,j)) / 2.0;
    const double vLeftDifference = (v(i-1,j) - v(i,j)) / 2.0;

    return (1.0/dx()) * (uRightMean*vRightMean - uLeftMean*vLeftMean) 
                + (alpha_/dx()) * (std::abs(uRightMean)*vRightDifference - std::abs(uLeftMean)*vLeftDifference);
}

double DonorCell::computeDuvDy(int i, int j) const
{
    const double vUpperMean = (v(i,j) + v(i+1,j)) / 2.0;
    const double uUpperMean = (u(i,j) + u(i,j+1)) / 2.0;
    const double vLowerMean = (v(i,j-1) + v(i+1,j-1)) / 2.0;
    const double uLowerMean = (u(i,j-1) + u(i,j)) / 2.0;

    const double uUpperDifference = (u(i,j) - u(i,j+1)) / 2.0;
    const double uLowerDifference = (u(i,j-1) - u(i,j)) / 2.0;
    

    return (1.0/dy()) * (vUpperMean*uUpperMean - vLowerMean*uLowerMean) 
            + (alpha_/dy()) * (std::abs(vUpperMean)*uUpperDifference - std::abs(vLowerMean)*uLowerDifference);
}