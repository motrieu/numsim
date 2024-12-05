#pragma once

#include <memory>
#include "pressureSolverParallel.h"

class SORParallel : public PressureSolverParallel
{

public:

    /// @brief constructor of SOR
    /// @param discretization shared pointer to the discretization which is either a central difference or donor cell object
    /// @param epsilon tolerance for the residuum for the iterative pressure solver
    /// @param maximumNumberOfIterations number of maximal iterations the pressure solver runs through
    /// @param partitioning contains information about the partition of the physical domain needed for MPI
    /// @param omega weight for correction term, for omega=1 one receives the Gauss-Seidel algorithm
    SORParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning, double omega);

    /// @brief implements the SOR-algorithm that solves for the new pressure p
    virtual void solve();

    virtual int getNumberOfIterations();

private:

    /// @brief computes pressure either on white or black tiles of the checker board used for parallelization of pressure communication
    /// @param leftAndLowerOffset switches between black and white tiles (if false the first tile is the one in the left-lower corner)
    void solveHalfStep(bool leftAndLowerOffset);

    /// @brief weight for correction term, for omega=1 one receives the Gauss-Seidel algorithm
    double omega_;

    int startPDifferences_;
    int numberIterations_;

};
