#include "pressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{
}

void PressureSolver::setBoundaryValues() 
{
    //sets boundary conditions for p: p(i,0) = p(i,1), p(i,N+1) = p(i,N)
    for (int i=(*discretization_).pIBegin()-1; i < (*discretization_).pIEnd()+1; i++)
    { 
        const double pInnerLower = (*discretization_).p(i,(*discretization_).pJBegin());
        const double pInnerUpper = (*discretization_).p(i,(*discretization_).pJEnd()-1);
        (*discretization_).p(i,(*discretization_).pJBegin()-1) = pInnerLower;
        (*discretization_).p(i,(*discretization_).pJEnd()) = pInnerUpper;
    }
    
    //sets boundary conditions for p: p(0,j) = p(1,j), p(N+1,j) = p(N,j)
    for (int j=(*discretization_).pJBegin()-1; j < (*discretization_).pJEnd()+1; j++)
    {
        const double pInnerLeft = (*discretization_).p((*discretization_).pIBegin(),j);
        const double pInnerRight = (*discretization_).p((*discretization_).pIEnd()-1,j);
        (*discretization_).p((*discretization_).pIBegin()-1,j) = pInnerLeft;
        (*discretization_).p((*discretization_).pIEnd(),j) = pInnerRight;
    }
}

const double PressureSolver::calcResNormSquared() const
{
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();
    
    double resNormSquared = 0;
    for (int i=(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
    {
        for (int j=(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
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