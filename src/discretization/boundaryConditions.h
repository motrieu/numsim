#pragma once

#include "staggeredGrid.h"
#include <array>


class BoundaryConditions : public StaggeredGrid
{
public:
    

    /// @brief constructor of boundary conditions
    /// @param nCells two-dimensional array for number of elements in x and y direction (halo cells not included)
    /// @param meshWidth two-dimensional array for mesh width in x and y direction
    BoundaryConditions(std::array<int,2> nCells, std::array<double,2> meshWidth);

    /// @brief perform NOSLIP condition for u and v in element i,j at the boundary indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 0 for north, 2 for east, 4 for south, 6 for west
    void noSlip(int i, int j, int directionIndex);

    /// @brief perform NOSLIP condition for u and v in element i,j at the corner indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 1 for north-east-corner, 3 for east-south-corner, 5 for south-west-corner, 7 for west-north-corner
    void noSlipCorner(int i, int j, int directionIndex);

    /// @brief perform NOSLIP condition for u and v in the diagonal neighbour of element i,j indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 1 for north-east-diagonal, 3 for east-south-diagonal, 5 for south-west-diagonal, 7 for west-north-diagonal 
    void noSlipDiagonal(int i, int j, int directionIndex);

    /// @brief perform SLIP condition for u and v in element i,j at the boundary indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 0 for north, 2 for east, 4 for south, 6 for west
    void slip(int i, int j, int directionIndex);

    /// @brief perform SLIP condition for u and v in the diagonal neighbour of element i,j indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 1 for north-east-diagonal, 3 for east-south-diagonal, 5 for south-west-diagonal, 7 for west-north-diagonal
    void slipDiagonal(int i, int j, int directionIndex);

    /// @brief perform INFLOW condition for u and v with predefined Dirichlet values uIN and vIN in element i,j at the boundary indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 0 for north, 2 for east, 4 for south, 6 for west 
    /// @param uIn Dirichlet value for u in element i,j
    /// @param vIn Dirichlet value for v in element i,j
    void inflow(int i, int j, int directionIndex, double uIn, double vIn);
    
    /// @brief perform INFLOW condition for u and v with Dirichlet values uIN and vIN in the diagonal neighbour of element i,j indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 1 for north-east-diagonal, 3 for east-south-diagonal, 5 for south-west-diagonal, 7 for west-north-diagonal
    /// @param uIn Dirichlet value for u
    /// @param vIn Dirichlet value for v
    void inflowDiagonal(int i, int j, int directionIndex, double uIn, double vIn);

    /// @brief perform OUTFLOW condition for u and v in element i,j at the boundary indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 0 for north, 2 for east, 4 for south, 6 for west
    void outflow(int i, int j, int directionIndex);

    /// @brief perform OUTFLOW condition for u and v in the diagonal neighbour of element i,j indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 1 for north-east-diagonal, 3 for east-south-diagonal, 5 for south-west-diagonal, 7 for west-north-diagonal
    void outflowDiagonal(int i, int j, int directionIndex);

    /// @brief perform PRESSURE condition for p with predefined Dirichlet value pRB  
    ///        either at the boundary of element i,j or in the diagonal neighbour indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 0 for north, 1 for north-east-diagonal, 2 for east, 3 for east-south-diagonal, 4 for south, 5 for south-west-diagonal, 6 for west, 7 for west-north-diagonal
    /// @param pRB Dirichlet value for p in element i,j
    void pressureDirichlet(int i, int j, int directionIndex, double pRB);

    /// @brief perform Neumann zero condition for p
    ///        either at the boundary of element i,j or in the diagonal neighbour indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 0 for north, 1 for north-east-diagonal, 2 for east, 3 for east-south-diagonal, 4 for south, 5 for south-west-diagonal, 6 for west, 7 for west-north-diagonal
    void pressureNeumannZero(int i, int j, int directionIndex);

    /// @brief perform Neumann zero condition for p in element i,j at the corner indicated with directionIndex
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @param directionIndex 1 for north-east-corner, 3 for east-south-corner, 5 for south-west-corner, 7 for west-north-corner
    void pressureNeumannZeroCorner(int i, int j, int directionIndex);
};