#include "pressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{
}

void PressureSolver::setBoundaryValues() 
{
    for (int i=0; i < (*discretization_).nCells()[0]; i++) //(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
    { 
        const double pInnerLower = (*discretization_).p(i,(*discretization_).pJBegin());
        const double pInnerUpper = (*discretization_).p(i,(*discretization_).pJEnd()-1);
        (*discretization_).p(i,(*discretization_).pJBegin()-1) = pInnerLower;
        (*discretization_).p(i,(*discretization_).pJEnd()) = pInnerUpper;
    }

    for (int j=0; j < (*discretization_).nCells()[1]; j++) //(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
    {
        const double pInnerLeft = (*discretization_).p((*discretization_).pIBegin(),j);
        const double pInnerRight = (*discretization_).p((*discretization_).pIEnd()-1,j);
        (*discretization_).p((*discretization_).pIBegin()-1,j) = pInnerLeft;
        (*discretization_).p((*discretization_).pIEnd(),j) = pInnerRight;
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