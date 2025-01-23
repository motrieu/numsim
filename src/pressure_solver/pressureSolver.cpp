#include "pressureSolver.h"


PressureSolver::PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations) :
    discretization_(discretization), epsilon_(epsilon), maximumNumberOfIterations_(maximumNumberOfIterations)
{
}

void PressureSolver::applyBoundaryConditions()
{
    for (int i = (*discretization_).setupIBegin(); i < (*discretization_).setupIEnd(); i++)
    {
        for (int j = (*discretization_).setupJBegin(); j < (*discretization_).setupJEnd(); j++)
        {
            if (((*discretization_).setup(i,j) != (*discretization_).indexFluid())
                    && (*discretization_).edgeDirections(i,j) != -1)
            {
                int edgeDirection = (*discretization_).edgeDirections(i,j);
                int numberFaces = (*discretization_).numberFaces(i,j);
                if ((edgeDirection%2 == 1) && (numberFaces == 2))
                {
                    (*discretization_).pressureNeumannZeroCorner(i, j, edgeDirection);
                }
                else
                {
                    if ((*discretization_).setup(i,j) == (*discretization_).indexNoSlip())
                        (*discretization_).pressureNeumannZero(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexSlip())
                        (*discretization_).pressureNeumannZero(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexInflow())
                        (*discretization_).pressureNeumannZero(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexOutflow())
                        (*discretization_).pressureNeumannZero(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexPressure())
                        (*discretization_).pressureDirichlet(i, j, edgeDirection, (*discretization_).pRB(i,j));
                }
            }
        }
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
            if ((*discretization_).setup(i,j) == (*discretization_).indexFluid())
            {
                const double rhs = (*discretization_).rhs(i,j);
                const double Pxx = ((*discretization_).p(i+1,j) - 2.0*(*discretization_).p(i,j) + (*discretization_).p(i-1,j)) / (dx*dx);
                const double Pyy = ((*discretization_).p(i,j+1) - 2.0*(*discretization_).p(i,j) + (*discretization_).p(i,j-1)) / (dy*dy);
                const double res = rhs - (Pxx + Pyy);

                resNormSquared += res*res;
            }
        }
    }

    return resNormSquared;
}