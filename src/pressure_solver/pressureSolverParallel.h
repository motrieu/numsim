#pragma once

#include <array>
#include <memory>
#include "pressureSolver.h"
#include "partitioning/partitioning.h"

class PressureSolverParallel : public PressureSolver
{
public:
    /// @brief constructor of pressure solver parallel
    /// @param discretization shared pointer to the discretization which is either a central difference or donor cell object
    /// @param epsilon tolerance for the residuum for the iterative pressure solver
    /// @param maximumNumberOfIterations number of maximal iterations the pressure solver runs through
    /// @param partitioning contains information about the partition of the physical domain needed for MPI
    PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning);

protected:

    /// @brief 
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