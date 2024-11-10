#pragma once

#include <memory>
#include "pressureSolver.h"

class SOR : public PressureSolver
{

public:

    /// @brief constructor of SOR
    /// @param discretization shared pointer to the discretization which is either a central difference or donor cell object
    /// @param epsilon tolerance for the residuum for the iterative pressure solver
    /// @param maximumNumberOfIterations number of maximal iterations the pressure solver runs through
    /// @param omega weight for correction term
    SOR(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, double omega);

    /// @brief implements the SOR-algorithm that solves for the new pressure p
    virtual void solve();

private:
    double omega_;

};
