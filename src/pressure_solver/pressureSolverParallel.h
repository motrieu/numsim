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

    /// @brief if the partition contains Dirichlet boundaries, then for each of them the boundary values for
    ///        the field variable p are set to the pressure values of the closest inner cells
    void setBoundaryValuesOnDirichletParallel();

    /// @brief sends the values needed by neigbouring processes after performing one half step of the pressure solver
    ///        receives the needed values for the next half step/iteration/time step from neighboring processes
    /// @param leftAndLowerOffset switches between black and white tiles (if false the first tile for the left and lower boundary is the one in the left-lower corner)
    /// @param rightOffset switches between black and white tiles (if false the first tile for the right boundary is the one in the right-lower corner)
    /// @param upperOffset switches between black and white tiles (if false the first tile for the upper boundary is the one in the left-upper corner)
    void receiveAndSendPressuresFromAndToOtherProcesses(bool leftAndLowerOffset, bool rightOffset, bool upperOffset);
    
    /// @brief calculates squared residual norm of pressure p, is used as termination criterium
    /// @return squared residual norm of pressure p reduced from all processes
    const double calcResNormSquaredParallel() const;


    /// @brief contains information about the partition of the physical domain needed for MPI
    Partitioning partitioning_;

    /// @brief defines whether left-lowest tile is black or white in the first half step (if false the first tile for the left and lower boundary is the one in the left-lower corner)
    bool leftAndLowerOffset_;

    /// @brief defines whether right-lowest tile is black or white in the first half step (if false the first tile for the right boundary is the one in the right-lower corner)
    bool rightOffset_;
    
    /// @brief defines whether left-upper tile is black or white in the first half step (if false the first tile for the upper boundary is the one in the left-upper corner)
    bool upperOffset_;

    const int pIBegin_;
    const int pIEnd_;
    const int pJBegin_;
    const int pJEnd_;

    const int nCellsX_;
    const int nCellsY_;
};