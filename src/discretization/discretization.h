#pragma once

#include "staggeredGrid.h"
#include <array>


class Discretization : public StaggeredGrid
{
public:
    
    /// @brief constructor of discretization
    /// @param nCells two-dimensional array for number of elements in x and y direction
    /// @param meshWidth two-dimensional array for mesh width in x and y direction
    Discretization(std::array<int,2> nCells, std::array<double,2> meshWidth);
  
    /// @brief compute the first derivative of u^2 in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(u^2)/dx in element i,j
    virtual double computeDu2Dx(int i, int j) const = 0;
 	
    /// @brief compute the first derivative of v^2 in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(v^2)/dy in element i,j
    virtual double computeDv2Dy(int i, int j) const = 0;
 	
    /// @brief compute the first derivative of u*v in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(uv)/dx in element i,j
    virtual double computeDuvDx(int i, int j) const = 0;
 	
    /// @brief compute the first derivative of u*v in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d(uv)/dy in element i,j
    virtual double computeDuvDy(int i, int j) const = 0;
 	
    /// @brief compute the second derivative of u in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d^2(u)/dx^2 in element i,j
    double computeD2uDx2(int i, int j) const;
 	
    /// @brief compute the second derivative of u in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d^2(u)/dy^2 in element i,j
    double computeD2uDy2(int i, int j) const;
 	
    /// @brief compute the second derivative of v in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d^2(v)/dx^2 in element i,j
    double computeD2vDx2(int i, int j) const;
 	
    /// @brief compute the second derivative of v in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return d^2(v)/dy^2 in element i,j
    double computeD2vDy2(int i, int j) const;
 	
    /// @brief compute the first derivative of p in x-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return dp/dx in element i,j
    double computeDpDx(int i, int j) const;
 	
    /// @brief compute the first derivative of p in y-direction in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return dp/dy in element i,j
    double computeDpDy(int i, int j) const;

};