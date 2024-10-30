#pragma once

#include <memory>
#include "pressure_solver/pressureSolver.h"

class GaussSeidel : public PressureSolver
{
    
public:

    /// @brief implements the Gauss-Seidel-algorithm that solves for the new pressure p
    virtual void solve();

};
