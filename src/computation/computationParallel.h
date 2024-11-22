// Inherits from computation
// Overloads runSimulation + additional protected methods
#pragma once

#include "computation.h"


class ComputationParallel : public Computation
{
public:
    void runSimulation();
    void initialize();

private:
    void applyBCOnDirichletBoundary();
    void applyPreliminaryBCOnDirichletBoundary();
    void applyBCInHaloCellsAtDirichletBoundary();
    void applyBCInHaloCellsAtNonDirichletBoundary();
    void applyBCOnNonDirichletBoundary();
    void applyPreliminaryBCOnNonDirichletBoundary();
    void computeTimeStepWidth();

    Partitioning partitioning_;
}