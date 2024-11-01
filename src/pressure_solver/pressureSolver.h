#pragma once

#include <memory>
#include "storage/discretization.h"

class PressureSolver
{
    
public:

    /// @brief constructor of pressure solver
    /// @param discretization shared pointer to the discretization which is either a central difference or donor cell object
    /// @param epsilon tolerance for the residuum for the iterative pressure solver
    /// @param maximumNumberOfIterations number of maximal iterations the pressure solver runs through
    PressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations);

    /// @brief solves for new pressure values, the new values are stored in the field variable p, implementation depends on which solver is used (GS or SOR)
    virtual void solve() = 0;

protected:

    /// @brief sets boundary values for field variable p, the pressure value of the closest inner cell is used to set the boundary 
    void setBoundaryValues();

    std::shared_ptr<Discretization> discretization_;
    double epsilon_;
    int	maximumNumberOfIterations_;

};
