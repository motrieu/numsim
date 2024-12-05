#pragma once

#include <array>
#include <memory>
#include <algorithm>

#include "settings/settings.h"
#include "discretization/discretization.h"
#include "discretization/donorCell.h"
#include "discretization/centralDifferences.h"
#include "pressure_solver/pressureSolver.h"
#include "pressure_solver/sor.h"
#include "pressure_solver/gaussSeidel.h"
#include "output_writer/output_writer_paraview.h"
#include "output_writer/output_writer_text.h"


class Computation
{
public:

    /// @brief initializes all objects needed in order to solve the Navier-Stokes equations
    ///        reads in the parameters from the given parameter file
    ///        decides based on parameter file which discretization scheme and which pressure solver to choose
    ///        applies the constant boundary conditions of u, v, F and G on the boundary faces
    /// @param argc number of arguments parsed in the command line, must be equal to 2
    /// @param argv the arguments parsed in the command line, contains the file name of the parameter file
    void initialize(int argc, char *argv[]);
 
    
    /// @brief solves the Navier-Stokes equations by looping over time and calculating for each time step the updated velocities u and v,  
    ///        pressure p, preliminary velocities F and G and the right hand side of the pressure Poisson problem (rhs)
    void runSimulation();

protected:

    /// @brief applies boundary conditions of u and v on the actual boundary faces
    ///        manipulates the in the data structure stored values such that they are overwritten
    ///        is only called once in the beginning of the computation
    ///        needs to be called before applyPreliminaryBCOnBoundary(), since it depends on the values set in applyBCOnBoundary()
    void applyBCOnBoundary();

    /// @brief applies boundary conditions of u and v in halo cells
    ///        manipulates the in the data structure stored values such that they are overwritten
    ///        has to be called in the beginning of every time step
    void applyBCInHaloCells();

    /// @brief applies boundary conditions of F and G on the actual boundary faces
    ///        manipulates the in the data structure stored values such that they are overwritten
    ///        is only called once in the beginning of the computation
    ///        the function applyBCOnBoundary() needs to be called prior to applyPreliminaryBCOnBoundary()
    void applyPreliminaryBCOnBoundary();
    
    /// @brief computes time step width for each time step such that the stability (diffusive and convective) is ensured
    void computeTimeStepWidth();
 	
    /// @brief computes F and G with possible external forces g[0], g[1] such as the gravitational force
    ///        manipulates the in the data structure stored values such that they are overwritten
    void computePreliminaryVelocities();
 	
    /// @brief computes right hand side (rhs) for the Poisson pressure equation
    ///        manipulates the in the data structure stored values such that they are overwritten
    void computeRightHandSide();

    /// @brief computes pressure by using the selected pressure solver algorithm (either Gauss-Seidel or SOR)
    ///        manipulates the in the data structure stored values such that they are overwritten
    void computePressure();
 	
    /// @brief computes new velocities u and v
    ///        manipulates the in the data structure stored values such that they are overwritten
    void computeVelocities();


    /// @brief contains parameters used for the current scenario
    Settings settings_;
 
    /// @brief shared pointer to the discretization, can either point towards the Central Difference or Donor Cell scheme
    std::shared_ptr<Discretization> discretization_;
    
    /// @brief shared pointer to the pressure solver, can either point towards the Gauss-Seidel or SOR algorithm
    std::unique_ptr<PressureSolver> pressureSolver_;
    
    /// @brief unique pointer to the output writer which produces vtk-output
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    
    /// @brief unique pointer to the output writer which produces readable txt-files
    std::unique_ptr<OutputWriterText> outputWriterText_;

    
    /// @brief array of doubles that contains the mesh width in x- and y-direction
    std::array<double,2> meshWidth_;
    
    /// @brief current time step width, is set in each time iteration
    double dt_;

    /// @brief mesh width in x-direction
    double dx_;

    /// @brief squared mesh width in x-direction
    double dxSquared_;

    /// @brief mesh width in y-direction
    double dy_;

    /// @brief squared mesh width in y-direction
    double dySquared_;
};
