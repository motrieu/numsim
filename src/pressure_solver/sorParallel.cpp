#include "sorParallel.h"

SORParallel::SORParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning, double omega) :
    PressureSolverParallel(discretization, epsilon, maximumNumberOfIterations, partitioning), omega_(omega)
{
}

void SORParallel::solve()
{
    setBoundaryValuesOnDirichletParallel();
    receiveAndSendPressuresFromAndToOtherProcesses();
}
