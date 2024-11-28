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
    /// @param secondHalfStep is false for the first half step of each iteration and true for the second
    void receiveAndSendPressuresFromAndToOtherProcesses(bool secondHalfStep);
    
    /// @brief calculates squared residual norm of pressure p, is used as termination criterium
    /// @return squared residual norm of pressure p reduced from all processes
    const double calcResNormSquaredParallel() const;


    /// @brief contains information about the partition of the physical domain needed for MPI
    Partitioning partitioning_;


    /// @brief defines whether left-lowest tile is black or white in the first half step (if false the first tile for the left and lower boundary is the one in the left-lower corner)
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<bool,2> leftAndLowerOffset_;

    /// @brief defines whether right-lowest tile is black or white in the first half step (if false the first tile for the right boundary is the one in the right-lower corner)
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> rightOffset_;

    /// @brief defines whether left-upper tile is black or white in the first half step (if false the first tile for the upper boundary is the one in the left-upper corner)
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> upperOffset_;


    /// @brief number of pressure values that need to be send to the left neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> sendLeftBufferLength_;

    /// @brief number of pressure values that need to be send to the lower neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> sendLowerBufferLength_;

    /// @brief number of pressure values that need to be send to the right neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> sendRightBufferLength_;

    /// @brief number of pressure values that need to be send to the upper neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> sendUpperBufferLength_;


    /// @brief number of pressure values that need to be received from the left neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> receiveLeftBufferLength_;

    /// @brief number of pressure values that need to be received from the lower neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> receiveLowerBufferLength_;

    /// @brief number of pressure values that need to be received from the right neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> receiveRightBufferLength_;

    /// @brief number of pressure values that need to be received from the upper neighbor after each half step of each iteration
    ///        first entry is valid for the first half step, the second entry for the second half step of each iteration
    std::array<int,2> receiveUpperBufferLength_;


    /// @brief first inner index for p in x-direction of the current process
    const int pIBegin_;

    /// @brief one after last inner index for p in x-direction of current process
    const int pIEnd_;

    /// @brief first inner index for p in y-direction of the current process
    const int pJBegin_;

    /// @brief one after last inner index for p in x-direction of current process
    const int pJEnd_;


    /// @brief number of elements in x direction for this process (halo cells not included)
    int nCellsX_;

    /// @brief number of elements in y direction for this process (halo cells not included)
    int nCellsY_;

};