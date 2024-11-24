#pragma once

#include <array>
#include <memory>
#include "pressureSolver.h"
#include "partitioning.h"

class PressureSolverParallel : public PressureSolver
{
public:
    PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning);

private:
    void setBoundaryValuesOnDirichletParallel();
    void receiveAndSendPressuresFromAndToOtherProcesses();
    const double calcResNormSquaredParallel() const;

    Partitioning partitioning_;

protected:
    bool startOffset_;
}