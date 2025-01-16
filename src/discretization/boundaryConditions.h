#pragma once

#include "staggeredGrid.h"
#include <array>


class BoundaryConditions : public StaggeredGrid
{
public:

    /// @brief 
    /// @param i 
    /// @param j 
    /// @param directionIndex 0 for north, 1 for east, 2 for south, 3 for west
    void noSlip(int i, int j, int directionIndex);
    /// @brief 
    /// @param i 
    /// @param j 
    /// @param directionIndex 0 for north-east-corner, 1 for east-south-corner, 2 for south-west-corner, 3 for west-north-corner
    void noSlipCorner(int i, int j, int directionIndex);
    void slip(int i, int j, int directionIndex);
    void inflow(int i, int j, int directionIndex, double uIn, double vIn);
    void outflow(int i, int j, int directionIndex);
    void pressureDirichlet(int i, int j, int directionIndex, double pRB);
    void pressureNeumannZero(int i, int j, int directionIndex);
    void pressureNeumannZeroCorner(int i, int j, int directionIndex);
};