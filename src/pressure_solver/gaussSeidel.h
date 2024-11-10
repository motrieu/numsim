#pragma once

#include "pressureSolver.h"

class GaussSeidel : public PressureSolver
{
    
public:

    using PressureSolver::PressureSolver;

    /// @brief implements the Gauss-Seidel-algorithm that solves for the new pressure p
    virtual void solve();

};
