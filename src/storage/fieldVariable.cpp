#include "fieldVariable.h"

#include <cmath>

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
  Array2D(size), origin_(origin), meshWidth_(meshWidth)
{
}

double FieldVariable::interpolateAt(double x, double y) const 
{
  const double xTransformed = (x - origin_[0])/meshWidth_[0];
  const double yTransformed = (y - origin_[1])/meshWidth_[1];

  const int leftXIndex = xTransformed;
  const int rightXIndex = leftXIndex + 1;
  const int lowerYIndex = yTransformed;
  const int upperYIndex = lowerYIndex + 1;

  const double percentageX = xTransformed - std::floor(xTransformed);
  const double percentageY = yTransformed - std::floor(yTransformed);

  const double fvLowerLeft = (*this)(leftXIndex, lowerYIndex);
  const double fvLowerRight = (*this)(rightXIndex, lowerYIndex);
  const double fvUpperLeft = (*this)(leftXIndex, upperYIndex);
  const double fvUpperRight = (*this)(rightXIndex, upperYIndex);
    
  const double horizontalInterpolationLower = fvLowerLeft + percentageX * (fvLowerRight - fvLowerLeft);
  const double horizontalInterpolationUpper = fvUpperLeft + percentageX * (fvUpperRight - fvUpperLeft);
  const double interpolatedValue = horizontalInterpolationLower + percentageY * (horizontalInterpolationUpper - horizontalInterpolationLower);

  return interpolatedValue;
}