#pragma once

#include "pressureSolver.h"

class GaussSeidel : public PressureSolver
{
    
public:

    // make PressureSolver constructor available in Gauss-Seidel
    using PressureSolver::PressureSolver;

    /// @brief implements the Gauss-Seidel-algorithm that solves for the new pressure p
    virtual void solve();

};
