#pragma once

#include "array2d.h"
#include <array>

/** A field variable is the discretization of a scalar function f(x) with x in the computational domain.
 *  More specifically, a scalar value is stored at discrete nodes/points. The nodes are arranged in an equidistant mesh
 *  with specified mesh width.
 */
class FieldVariable : public Array2D
{
public:
    
    /// @brief constructor of field variable
    /// @param size two-dimensional array for number of elements in x and y direction (halo cells included)
    /// @param origin cartesian coordinates of the point with (i,j) = (0,0), this is different from (0,0) for all field variables 
    /// @param meshWidth two-dimensional array for mesh width in x and y direction
    FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth);

    /// @brief get the value of the field variable at the cartesian coordinate (x,y)
    ///        the value is bilinearly interpolated between adjacent field variable values
    /// @param x cartesian x-coordinate
    /// @param y cartesian y-coordinate
    /// @return inerpolated value at cartesian coordinate (x,y)
    double interpolateAt(double x, double y) const;

private:
    /// @brief cartesian coordinates of the point with (i,j) = (0,0), this is different from (0,0) for all field variables 
    const std::array<double,2> origin_;

    /// @brief array of doubles that contains the mesh width in x- and y-direction
    const std::array<double,2> meshWidth_;
};