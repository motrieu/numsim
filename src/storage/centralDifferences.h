#pragma once

#include "discretization.h"

class CentralDifferences : public Discretization
{
public:
    using Discretization::Discretization;

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
};