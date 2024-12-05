#include "gaussSeidel.h"

void GaussSeidel::solve()
{
    int n = 0;
    double resNormSquared;

    do
    {
        for (int i=(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
        {
            for (int j=(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
            {
                const double prefactor = (dxSquared_ * dySquared_) / (2.0*(dxSquared_ + dySquared_));
                const double firstSummand = ((*discretization_).p(i-1,j) + (*discretization_).p(i+1,j)) / dxSquared_;
                const double secondSummand = ((*discretization_).p(i,j-1) + (*discretization_).p(i,j+1)) / dySquared_;
                const double rhs = (*discretization_).rhs(i,j);

                (*discretization_).p(i,j) = prefactor * (firstSummand + secondSummand - rhs);
            }
        }
        
        resNormSquared = calcResNormSquared();
        n++;

        setBoundaryValues();
    }
    //Termination criteria: either number of maximal iterations is reached or residual squared norm is less or equal to given threshold
    while ((n < maximumNumberOfIterations_) && (resNormSquared > numberOfValues_*epsilonSquared_));
}