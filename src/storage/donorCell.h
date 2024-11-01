#pragma once

#include "discretization.h"


class DonorCell : public Discretization
{
public:

    /// @brief constructor of donor cell discretization scheme
    /// @param nCells two-dimensional array for number of elements in x and y direction
    /// @param meshWidth two-dimensional array for mesh width in x and y direction
    /// @param alpha alpha-parameter which maps between central differences and donor cell, 0=<alpha=<1, for alpha=0 the scheme is the pure central differences scheme, for alpha=1 the pure donor cell scheme
    DonorCell(std::array<int,2> nCells, std::array<double,2> meshWidth, double alpha);
    
    /// @brief compute the first derivative of u^2 in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(u^2)/dx in element i,j
    virtual double computeDu2Dx(int i, int j) const;
 	
    /// @brief compute the first derivative of v^2 in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(v^2)/dy in element i,j
    virtual double computeDv2Dy(int i, int j) const;
 	
    /// @brief compute the first derivative of u*v in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(uv)/dx in element i,j
    virtual double computeDuvDx(int i, int j) const;
 	
    /// @brief compute the first derivative of u*v in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(uv)/dy in element i,j
    virtual double computeDuvDy(int i, int j) const;

private:
    double alpha_;

};