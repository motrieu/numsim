#include "sorParallel.h"

#include <mpi.h>

SORParallel::SORParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning, double omega) :
    PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning), omega_(omega)
{
}

void SORParallel::solve()
{
    int n = 0;
    double resNormSquaredParallel;
    const double epsSquared = epsilon_*epsilon_;
    const double numberOfValues = ((*discretization_).nCells()[0]) * ((*discretization_).nCells()[1]);

    do
    {
        solveHalfStep(leftAndLowerOffset_);
        receiveAndSendPressuresFromAndToOtherProcesses(leftAndLowerOffset_, rightOffset_, upperOffset_);
        solveHalfStep(!leftAndLowerOffset_);
        receiveAndSendPressuresFromAndToOtherProcesses(!leftAndLowerOffset_, !rightOffset_, !upperOffset_);

        resNormSquaredParallel = calcResNormSquaredParallel();
        n++;

        setBoundaryValuesOnDirichletParallel();
    }
    //Termination criteria: either number of maximal iterations is reached or residual squared norm is less or equal to given threshold
    while ((n < maximumNumberOfIterations_) && (resNormSquaredParallel > numberOfValues*epsSquared));
}

void SORParallel::solveHalfStep(bool leftAndLowerOffset)
{
    const double dx = (*discretization_).dx();
    const double dy = (*discretization_).dy();

    bool rowOffset = leftAndLowerOffset;
    for (int j=pJBegin_; j < pJEnd_; j++)
    {
        for (int i=pIBegin_+rowOffset; i < pIEnd_; i+=2)
        {
            const double prefactor = (dx*dx*dy*dy) / (2.0*(dx*dx+dy*dy));
            const double firstSummand = ((*discretization_).p(i-1,j) + (*discretization_).p(i+1,j)) / (dx*dx);
            const double secondSummand = ((*discretization_).p(i,j-1) + (*discretization_).p(i,j+1)) / (dy*dy);
            const double rhs = (*discretization_).rhs(i,j);
            const double p = (*discretization_).p(i,j);

            (*discretization_).p(i,j) = p + omega_ * (prefactor * (firstSummand + secondSummand - rhs) - p);
        }
        rowOffset = !rowOffset;
    }
}