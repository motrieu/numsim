// Inherits from computation
// Overloads runSimulation + additional protected methods
#pragma once

#include "computation.h"


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
    

    

    Partitioning partitioning_;
};