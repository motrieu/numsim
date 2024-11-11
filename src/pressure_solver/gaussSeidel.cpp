#include "gaussSeidel.h"

void GaussSeidel::solve()
{
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();

    int n = 0;
    double resNormSquared;
    const double epsSquared = epsilon_*epsilon_;
    const double numberOfValues = ((*discretization_).nCells()[0]-2) * ((*discretization_).nCells()[1]-2);
    
    do
    {
        for (int i=(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
        {
            for (int j=(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
            {
                const double prefactor = (dx*dx*dy*dy) / (2.0*(dx*dx+dy*dy));
                const double firstSummand = ((*discretization_).p(i-1,j) + (*discretization_).p(i+1,j)) / (dx*dx);
                const double secondSummand = ((*discretization_).p(i,j-1) + (*discretization_).p(i,j+1)) / (dy*dy);
                const double rhs = (*discretization_).rhs(i,j);

                (*discretization_).p(i,j) = prefactor * (firstSummand + secondSummand - rhs);
            }
        }
        
        resNormSquared = calcResNormSquared();
        n++;

        setBoundaryValues();
    }
    //Termination criteria: either number of maximal iterations is reached or residual squared norm is less or equal to given threshold
    while ((n < maximumNumberOfIterations_) && (resNormSquared > numberOfValues*epsSquared));
}