#include "boundaryConditions.h"

BoundaryConditions::BoundaryConditions(std::array<int,2> nCells, std::array<double,2> meshWidth) :
    StaggeredGrid(nCells, meshWidth)
{
}

void BoundaryConditions::noSlip(int i, int j, int directionIndex)
{
    if (directionIndex == 0)
    {
        v(i,j) = 0.0;
        u(i,j) = -u(i,j+1);
    }
    else if (directionIndex == 2)
    {
        u(i,j) = 0.0;
        v(i,j) = -v(i+1,j);
    }
    else if (directionIndex == 4)
    {
        v(i,j-1) = 0.0;
        u(i,j) = -u(i,j-1);
    }
    else if (directionIndex == 6)
    {
        u(i-1,j) = 0.0;
        v(i,j) = -v(i-1,j);
    }
}

void BoundaryConditions::noSlipCorner(int i, int j, int directionIndex)
{
    if (directionIndex == 1)
    {
        u(i,j) = 0.0;
        v(i,j) = 0.0;   
    }
    else if (directionIndex == 3)
    {
        u(i,j) = 0.0;
        v(i,j-1) = 0.0;
    }
    else if (directionIndex == 5)
    {
        u(i-1,j) = 0.0;
        v(i,j-1) = 0.0;
    }
    else if (directionIndex == 7)
    {
        u(i-1,j) = 0.0;
        v(i,j) = 0.0;
    }
}

void BoundaryConditions::noSlipDiagonal(int i, int j, int directionIndex)
{
    if (directionIndex == 1)
    {
        u(i,j) = -u(i,j+1);
        v(i,j) = -v(i+1,j);
    }
    else if (directionIndex == 3)
    {
        u(i,j) = -u(i,j-1);
        v(i,j-1) = -v(i+1,j-1);
    }
    else if (directionIndex == 5)
    {
        u(i-1,j) = -u(i-1,j-1);
        v(i,j-1) = -v(i-1,j-1);
    }
    else if (directionIndex == 7)
    {
        u(i-1,j) = -u(i-1,j+1);
        v(i,j) = -v(i-1,j);
    }
}

void BoundaryConditions::slip(int i, int j, int directionIndex)
{
    if (directionIndex == 0)
    {
        v(i,j) = 0.0;
        u(i,j) = u(i,j+1);
    }
    else if (directionIndex == 2)
    {
        u(i,j) = 0.0;
        v(i,j) = v(i+1,j);
    }
    else if (directionIndex == 4)
    {
        v(i,j-1) = 0.0;
        u(i,j) = u(i,j-1);
    }
    else if (directionIndex == 6)
    {
        u(i-1,j) = 0.0;
        v(i,j) = v(i-1,j);
    }
}

void BoundaryConditions::slipDiagonal(int i, int j, int directionIndex)
{
    if (directionIndex == 1)
    {
        u(i,j) = u(i,j+1);
        v(i,j) = v(i+1,j);
    }
    else if (directionIndex == 3)
    {
        u(i,j) = u(i,j-1);
        v(i,j-1) = v(i+1,j-1);
    }
    else if (directionIndex == 5)
    {
        u(i-1,j) = u(i-1,j-1);
        v(i,j-1) = v(i-1,j-1);
    }
    else if (directionIndex == 7)
    {
        u(i-1,j) = u(i-1,j+1);
        v(i,j) = v(i-1,j);
    }
}

void BoundaryConditions::inflow(int i, int j, int directionIndex, double uIn, double vIn)
{
    if (directionIndex == 0)
    {
        v(i,j) = vIn;
        u(i,j) = 2*uIn - u(i,j+1);
    }
    else if (directionIndex == 2)
    {
        u(i,j) = uIn;
        v(i,j) = 2*vIn - v(i+1,j);
    }
    else if (directionIndex == 4)
    {
        v(i,j-1) = vIn;
        u(i,j) = 2*uIn - u(i,j-1);
    }
    else if (directionIndex == 6)
    {
        u(i-1,j) = uIn;
        v(i,j) = 2*vIn - v(i-1,j);
    }
}

void BoundaryConditions::inflowDiagonal(int i, int j, int directionIndex, double uIn, double vIn)
{
    if (directionIndex == 1)
    {
        u(i,j) = 2*uIn - u(i,j+1);
        v(i,j) = 2*vIn - v(i+1,j);
    }
    else if (directionIndex == 3)
    {
        u(i,j) = 2*uIn - u(i,j-1);
        v(i,j-1) = 2*vIn - v(i+1,j-1);
    }
    else if (directionIndex == 5)
    {
        u(i-1,j) = 2*uIn - u(i-1,j-1);
        v(i,j-1) = 2*vIn - v(i-1,j-1);
    }
    else if (directionIndex == 7)
    {
        u(i-1,j) = 2*uIn - u(i-1,j+1);
        v(i,j) = 2*vIn - v(i-1,j);
    }
}

void BoundaryConditions::outflow(int i, int j, int directionIndex)
{
    if (directionIndex == 0)
    {
        u(i,j) = u(i,j+1);
        v(i,j) = v(i,j+1);
    }
    else if (directionIndex == 2)
    {
        u(i,j) = u(i+1,j);
        v(i,j) = v(i+1,j);
    }
    else if (directionIndex == 4)
    {
        u(i,j) = u(i,j-1);
        v(i,j-1) = v(i,j-2);
    }
    else if (directionIndex == 6)
    {
        u(i-1,j) = u(i-2,j);
        v(i,j) = v(i-1,j);
    }
}

void BoundaryConditions::outflowDiagonal(int i, int j, int directionIndex)
{
    if (directionIndex == 1)
    {
        u(i,j) = u(i,j+1);
        v(i,j) = v(i+1,j);
    }
    else if (directionIndex == 3)
    {
        u(i,j) = u(i,j-1);
        v(i,j-1) = v(i+1,j-1);
    }
    else if (directionIndex == 5)
    {
        u(i-1,j) = u(i-1,j-1);
        v(i,j-1) = v(i-1,j-1);
    }
    else if (directionIndex == 7)
    {
        u(i-1,j) = u(i-1,j+1);
        v(i,j) = v(i-1,j);
    }
}

void BoundaryConditions::pressureDirichlet(int i, int j, int directionIndex, double pRB)
{
    if (directionIndex == 0)
    {
        p(i,j) = 2*pRB - p(i,j+1);
    }
    if (directionIndex == 1)
    {
        p(i,j) = 4*pRB - p(i,j+1) - p(i+1,j+1) - p(i+1,j);
    }
    else if (directionIndex == 2)
    {
        p(i,j) = 2*pRB - p(i+1,j);
    }
    else if (directionIndex == 3)
    {
        p(i,j) = 4*pRB - p(i+1,j) - p(i+1,j-1) - p(i,j-1);
    }
    else if (directionIndex == 4)
    {
        p(i,j) = 2*pRB - p(i,j-1);
    }
    else if (directionIndex == 5)
    {
        p(i,j) = 4*pRB - p(i,j-1) - p(i-1,j-1) - p(i-1,j);
    }
    else if (directionIndex == 6)
    {
        p(i,j) = 2*pRB - p(i-1,j);
    }
    else if (directionIndex == 7)
    {
        p(i,j) = 4*pRB - p(i-1,j) - p(i-1,j+1) - p(i,j+1);
    }
}

void BoundaryConditions::pressureNeumannZero(int i, int j, int directionIndex)
{
    if (directionIndex == 0)
    {
        p(i,j) = p(i,j+1);
    }
    else if (directionIndex == 1)
    {
        p(i,j) = 3*p(i+1,j+1) - p(i,j+1) - p(i+1,j);
    }
    else if (directionIndex == 2)
    {
        p(i,j) = p(i+1,j);
    }
    else if (directionIndex == 3)
    {
        p(i,j) = 3*p(i+1,j-1) - p(i+1,j) - p(i,j-1);
    }
    else if (directionIndex == 4)
    {
        p(i,j) = p(i,j-1);
    }
    else if (directionIndex == 5)
    {
        p(i,j) = 3*p(i-1,j-1) - p(i,j-1) - p(i-1,j);
    }
    else if (directionIndex == 6)
    {
        p(i,j) = p(i-1,j);
    }
    else if (directionIndex == 7)
    {
        p(i,j) = 3*p(i-1,j+1) - p(i-1,j) - p(i,j+1);
    }
}

void BoundaryConditions::pressureNeumannZeroCorner(int i, int j, int directionIndex)
{
    if (directionIndex == 1)
    {
        p(i,j) = 0.5 * (p(i,j+1) + p(i+1,j));
    }
    else if (directionIndex == 3)
    {
        p(i,j) = 0.5 * (p(i+1,j) + p(i,j-1));
    }
    else if (directionIndex == 5)
    {
        p(i,j) = 0.5 * (p(i,j-1) + p(i-1,j));
    }
    else if (directionIndex == 7)
    {
        p(i,j) = 0.5 * (p(i-1,j) + p(i,j+1));
    }
}
