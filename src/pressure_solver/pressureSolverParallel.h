#pragma once

#include <array>
#include <memory>
#include "pressureSolver.h"
#include "partitioning/partitioning.h"

class PressureSolverParallel : public PressureSolver
{
public:
    PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning);

protected:
    void setBoundaryValuesOnDirichletParallel();
    void receiveAndSendPressuresFromAndToOtherProcesses();
    const double calcResNormSquaredParallel() const;

    Partitioning partitioning_;
    bool startOffset_;
};