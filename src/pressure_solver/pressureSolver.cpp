#include "pressure_solver/pressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{
}

void PressureSolver::setBoundaryValues() 
{
    const int nCellsX = (*discretization_).nCells()[0];
    const int nCellsY = (*discretization_).nCells()[1];

    for (int i=1; i < nCellsX-1; i++)
    { 
        const double pInnerLower = (*discretization_).p(i,1);
        const double pInnerUpper = (*discretization_).p(i,nCellsY-2);
        (*discretization_).p(i,0) = pInnerLower;
        (*discretization_).p(i,nCellsY-1) = pInnerUpper;
    }

    for (int j=1; j < nCellsY-1; j++)
    {
        const double pInnerLeft = (*discretization_).p(1,j);
        const double pInnerRight = (*discretization_).p(nCellsX-2,j);
        (*discretization_).p(0,j) = pInnerLeft;
        (*discretization_).p(nCellsX-1,j) = pInnerRight;
    }
    
}