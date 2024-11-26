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
    void receiveAndSendPressuresFromAndToOtherProcesses(bool leftAndLowerOffset, bool rightOffset, bool upperOffset);
    const double calcResNormSquaredParallel() const;

    Partitioning partitioning_;

    bool leftAndLowerOffset_;
    bool rightOffset_;
    bool upperOffset_;

    const int pIBegin_;
    const int pIEnd_;
    const int pJBegin_;
    const int pJEnd_;

    const int nCellsX_;
    const int nCellsY_;
};