#include "pressure_solver/gaussSeidel.h"

void GaussSeidel::solve()
{
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();

    const int pIBegin = (*discretization_).pIBegin();
    const int pJBegin = (*discretization_).pJBegin();
    const int pIEnd = (*discretization_).pIEnd();
    const int pJEnd = (*discretization_).pJEnd();

    int n = 0;
    double resNormSquared = 1.0;
    const double epsSquared = epsilon_*epsilon_;
    const double numberOfValues = (*discretization_).nCells()[0] * (*discretization_).nCells()[1];
    
    while ((n < maximumNumberOfIterations_) && (resNormSquared > numberOfValues*epsSquared))
    {
        for (int j=pJBegin; j < pJEnd; j++)
        {
            for (int i=pIBegin; i < pIEnd; i++)
            {
                const double prefactor = (dx*dx*dy*dy) / (2.0*(dx*dx+dy*dy));
                const double firstSummand = ((*discretization_).p(i-1,j) + (*discretization_).p(i+1,j)) / (dx*dx);
                const double secondSummand = ((*discretization_).p(i,j-1) + (*discretization_).p(i,j+1)) / (dy*dy);
                const double rhs = (*discretization_).rhs(i,j);

                (*discretization_).p(i,j) = prefactor * (firstSummand + secondSummand - rhs);
            }
        }
        
        resNormSquared = 0;
         for (int j=pJBegin; j < pJEnd; j++)
        {
            for (int i=pIBegin; i < pIEnd; i++)
            {
                const double rhs = (*discretization_).rhs(i,j);
                const double Pxx = ((*discretization_).p(i+1,j) - 2.0*(*discretization_).p(i,j) + (*discretization_).p(i-1,j)) / (dx*dx);
                const double Pyy = ((*discretization_).p(i,j+1) - 2.0*(*discretization_).p(i,j) + (*discretization_).p(i,j-1)) / (dy*dy);
                const double res = rhs - (Pxx + Pyy);
                resNormSquared += res*res;
            }
        }
    }
}