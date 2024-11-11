#include "fieldVariable.h"

#include <cmath>

FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
  Array2D(size), origin_(origin), meshWidth_(meshWidth)
{
}

double FieldVariable::interpolateAt(double x, double y) const 
{
  // consider that different field variables live on different parts of the cell, which is reflected in the different origins
  const double xTransformed = (x - origin_[0])/meshWidth_[0];
  const double yTransformed = (y - origin_[1])/meshWidth_[1];

  // obtain the adjacent cell indices by integer casting 
  const int leftXIndex = xTransformed;
  const int rightXIndex = leftXIndex + 1;
  const int lowerYIndex = yTransformed;
  const int upperYIndex = lowerYIndex + 1;

  // calculates how to weight the interpolations in x- and y-direction 
  const double percentageX = xTransformed - std::floor(xTransformed);
  const double percentageY = yTransformed - std::floor(yTransformed);

  const double fvLowerLeft = (*this)(leftXIndex, lowerYIndex);
  const double fvLowerRight = (*this)(rightXIndex, lowerYIndex);
  const double fvUpperLeft = (*this)(leftXIndex, upperYIndex);
  const double fvUpperRight = (*this)(rightXIndex, upperYIndex);
  
  // first: linear interpolation in x-direction
  const double horizontalInterpolationLower = fvLowerLeft + percentageX * (fvLowerRight - fvLowerLeft);
  const double horizontalInterpolationUpper = fvUpperLeft + percentageX * (fvUpperRight - fvUpperLeft);

  // second: linear interpolation in y-direction between the interpolated x-values
  const double interpolatedValue = horizontalInterpolationLower + percentageY * (horizontalInterpolationUpper - horizontalInterpolationLower);

  return interpolatedValue;
}