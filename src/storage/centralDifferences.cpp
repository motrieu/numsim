#include "storage/centralDifferences.h"

double CentralDifferences::computeDu2Dx(int i, int j) const
{
    double rightMean = (u(i,j) + u(i+1,j)) / 2.0;
    double leftMean = (u(i-1,j) + u(i,j)) / 2.0;

    return (rightMean*rightMean - leftMean*leftMean) / dx();
}
 	
double CentralDifferences::computeDv2Dy(int i, int j) const
{
    double upperMean = (v(i,j) + v(i,j+1)) / 2.0;
    double lowerMean = (v(i,j-1) + v(i,j)) / 2.0;

    return (upperMean*upperMean - lowerMean*lowerMean) / dy();
}

double CentralDifferences::computeDuvDx(int i, int j) const
{
    double vRightMean = (v(i,j) + v(i+1,j)) / 2.0;
    double uRightMean = (u(i,j) + u(i,j+1)) / 2.0;
    double vLeftMean = (v(i-1,j) + v(i,j)) / 2.0;
    double uLeftMean = (u(i-1,j) + u(i-1,j+1)) / 2.0;

    return (vRightMean*uRightMean - vLeftMean*uLeftMean) / dx();
}

double CentralDifferences::computeDuvDy(int i, int j) const
{
    double uUpperMean = (u(i,j) + u(i,j+1)) / 2.0;
    double vUpperMean = (v(i,j) + v(i+1,j)) / 2.0;
    double uLowerMean = (u(i,j-1) + u(i,j)) / 2.0;
    double vLowerMean = (v(i,j-1) + v(i+1,j-1)) / 2.0;

    return (uUpperMean*vUpperMean - uLowerMean*vLowerMean) / dy();
}