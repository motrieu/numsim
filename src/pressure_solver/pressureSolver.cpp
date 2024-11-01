#include "pressureSolver.h"


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

const double PressureSolver::calc2NormOfP() const
{
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();
    
    double resNormSquared = 0;
    for (int j=(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
    {
        for (int i=(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
        {
            const double rhs = (*discretization_).rhs(i,j);
            const double Pxx = ((*discretization_).p(i+1,j) - 2.0*(*discretization_).p(i,j) + (*discretization_).p(i-1,j)) / (dx*dx);
            const double Pyy = ((*discretization_).p(i,j+1) - 2.0*(*discretization_).p(i,j) + (*discretization_).p(i,j-1)) / (dy*dy);
            const double res = rhs - (Pxx + Pyy);
            resNormSquared += res*res;
        }
    }
    return resNormSquared;
}