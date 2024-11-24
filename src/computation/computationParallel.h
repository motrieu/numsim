// Inherits from computation
// Overloads runSimulation + additional protected methods
#pragma once

#include "computation.h"
#include "pressure_solver/pressureSolverParallel.h"
#include "pressure_solver/sorParallel.h"

class ComputationParallel : public Computation
{
public:
    void runSimulation();
    void initialize(int argc, char *argv[]);

private:
    void applyBCOnDirichletBoundary();
    void applyPreliminaryBCOnDirichletBoundary();
    void applyBCInHaloCellsAtDirichletBoundary();
    void receiveAndSendVelocitiesFromAndToOtherProcesses();
    void computeTimeStepWidthParallel();
    void computePreliminaryVelocities();
    void receiveAndSendPreliminaryVelocitiesFromAndToOtherProcesses();
    

    std::unique_ptr<PressureSolverParallel> pressureSolverParallel_;

    Partitioning partitioning_;
};