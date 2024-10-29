#include "storage/fieldVariable.h"

#include <cmath>

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
  Array2D(size), origin_(origin), meshWidth_(meshWidth)
{
}

double FieldVariable::interpolateAt(double x, double y) const 
{
    double xTransformed = (x - origin_[0])/meshWidth_[0];
    double yTransformed = (y - origin_[1])/meshWidth_[1];

    int leftXIndex = xTransformed;
    int rightXIndex = leftXIndex + 1;
    int lowerYIndex = yTransformed;
    int upperYIndex = lowerYIndex + 1;

    double percentageX = xTransformed - std::floor(xTransformed);
    double percentageY = yTransformed - std::floor(yTransformed);

    double fvLowerLeft = (*this)(leftXIndex, lowerYIndex);
    double fvLowerRight = (*this)(rightXIndex, lowerYIndex);
    double fvUpperLeft = (*this)(leftXIndex, upperYIndex);
    double fvUpperRight = (*this)(rightXIndex, upperYIndex);
    
    double horizontalInterpolationLower = fvLowerLeft + percentageX * (fvLowerRight - fvLowerLeft);
    double horizontalInterpolationUpper = fvUpperLeft + percentageX * (fvUpperRight - fvUpperLeft);
    double interpolatedValue = horizontalInterpolationLower + percentageY * (horizontalInterpolationUpper - horizontalInterpolationLower);

    return interpolatedValue;
}