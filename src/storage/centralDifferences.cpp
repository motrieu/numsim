#include "centralDifferences.h"

double CentralDifferences::computeDu2Dx(int i, int j) const
{
    const double uRightMean = (u(i,j) + u(i+1,j)) / 2.0;
    const double uLeftMean = (u(i-1,j) + u(i,j)) / 2.0;

    return (uRightMean*uRightMean - uLeftMean*uLeftMean) / dx();
}
 	
double CentralDifferences::computeDv2Dy(int i, int j) const
{
    const double vUpperMean = (v(i,j) + v(i,j+1)) / 2.0;
    const double vLowerMean = (v(i,j-1) + v(i,j)) / 2.0;

    return (vUpperMean*vUpperMean - vLowerMean*vLowerMean) / dy();
}

double CentralDifferences::computeDuvDx(int i, int j) const
{
    const double vRightMean = (v(i,j) + v(i+1,j)) / 2.0;
    const double uRightMean = (u(i,j) + u(i,j+1)) / 2.0;
    const double vLeftMean = (v(i-1,j) + v(i,j)) / 2.0;
    const double uLeftMean = (u(i-1,j) + u(i-1,j+1)) / 2.0;

    return (vRightMean*uRightMean - vLeftMean*uLeftMean) / dx();
}

double CentralDifferences::computeDuvDy(int i, int j) const
{
    const double uUpperMean = (u(i,j) + u(i,j+1)) / 2.0;
    const double vUpperMean = (v(i,j) + v(i+1,j)) / 2.0;
    const double uLowerMean = (u(i,j-1) + u(i,j)) / 2.0;
    const double vLowerMean = (v(i,j-1) + v(i+1,j-1)) / 2.0;

    return (uUpperMean*vUpperMean - uLowerMean*vLowerMean) / dy();
}