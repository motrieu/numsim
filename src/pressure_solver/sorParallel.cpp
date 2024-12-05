#include "sorParallel.h"

#include <mpi.h>

SORParallel::SORParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning, double omega) :
    PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning), omega_(omega)
{
    startPDifferences_ = 2;
}

void SORParallel::solve()
{
    int n = 0;
    double resNormSquaredParallel;

    if (startPDifferences_ == 0)
    {
        for (int i=pIBegin_; i < pIEnd_; i++)
        {
            for (int j=pJBegin_; j < pJEnd_; j++)
            {
                const double p = (*discretization_).p(i,j);
                const double pOld = (*discretization_).pOld(i,j);
                const double pOld2 = (*discretization_).pOld2(i,j);
                (*discretization_).p(i,j) = p + 0.8 * (1.5*p - pOld + 0.5*pOld2);
                (*discretization_).pOld(i,j) = p;
                (*discretization_).pOld2(i,j) = pOld;
            }
        }
    }
    else
    {
        for (int i=pIBegin_; i < pIEnd_; i++)
        {
            for (int j=pJBegin_; j < pJEnd_; j++)
            {
                const double p = (*discretization_).p(i,j);
                const double pOld = (*discretization_).pOld(i,j);
                (*discretization_).pOld(i,j) = p;
                (*discretization_).pOld2(i,j) = pOld;
            }
        }
        startPDifferences_--;
    }

    do
    {
        // solve first half of the checker board (either white or black tiles depending on position of partition in global context)
        solveHalfStep(leftAndLowerOffset_[0]);

        // comunicate calculated pressures to neighboring processes
        receiveAndSendPressuresFromAndToOtherProcesses(false);

        // solve second half of the checker board (either black or whiles tiles (opposite to prior half step) depending on position of partition in global context)
        solveHalfStep(leftAndLowerOffset_[1]);

        // comunicate calculated pressures to neighboring processes such that in the next iteration/time step everything has been updated
        receiveAndSendPressuresFromAndToOtherProcesses(true);

        if (n%10 == 0)
            resNormSquaredParallel = calcResNormSquaredParallel();
        n++;

        setBoundaryValuesOnDirichletParallel();
    }
    //Termination criteria: either number of maximal iterations is reached or residual squared norm is less or equal to given threshold
    while ((n < maximumNumberOfIterations_) && (resNormSquaredParallel > numberOfValuesGlobal_*epsilonSquared_));
    
    numberIterations_ = n;
}

void SORParallel::solveHalfStep(bool leftAndLowerOffset)
{
    bool rowOffset = leftAndLowerOffset;
    for (int j=pJBegin_; j < pJEnd_; j++)
    {
        for (int i=pIBegin_+rowOffset; i < pIEnd_; i+=2)
        {
            const double prefactor = (dxSquared_*dySquared_) / (2.0*(dxSquared_+dySquared_));
            const double firstSummand = ((*discretization_).p(i-1,j) + (*discretization_).p(i+1,j)) / dxSquared_;
            const double secondSummand = ((*discretization_).p(i,j-1) + (*discretization_).p(i,j+1)) / dySquared_;
            const double rhs = (*discretization_).rhs(i,j);
            const double p = (*discretization_).p(i,j);

            (*discretization_).p(i,j) = p + omega_ * (prefactor * (firstSummand + secondSummand - rhs) - p);
        }
        rowOffset = !rowOffset;
    }
}

int SORParallel::getNumberOfIterations()
{
    return numberIterations_;
}