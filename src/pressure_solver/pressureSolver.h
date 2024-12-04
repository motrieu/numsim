#pragma once

#include <memory>
#include "discretization/discretization.h"

class PressureSolver
{
    
public:

    /// @brief constructor of pressure solver
    /// @param discretization shared pointer to the discretization which is either a central difference or donor cell object
    /// @param epsilon tolerance for the residuum for the iterative pressure solver
    /// @param maximumNumberOfIterations number of maximal iterations the pressure solver runs through
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /// @brief solves for new pressure values, the new values are stored in the field variable p, implementation depends on which solver is used (GS or SOR)
    virtual void solve(int resNormIntervall) = 0;

protected:

    /// @brief sets boundary values for field variable p, the pressure value of the closest inner cell is used to set the boundary 
    void setBoundaryValues();

    /// @brief calculates squared residual norm of pressure p, is used as termination criterium
    /// @return squared residual norm of pressure p in the current process
    const double calcResNormSquared() const;

    /// @brief shared pointer to the discretization, can either point towards the Central Difference or Donor Cell scheme
    std::shared_ptr<Discretization> discretization_;

    /// @brief squared threshold that is used for the termination criterium in order to find out if the solver has converged
    double epsilonSquared_;

    /// @brief maximal number of iterations used as termination criterium
    int	maximumNumberOfIterations_;

    /// @brief squared mesh width in x-direction
    double dxSquared_;

    /// @brief squared mesh width in y-direction
    double dySquared_;

    /// @brief number of cells of the physical domain of the current process
    double numberOfValues_;

};
