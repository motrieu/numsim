// Inherits from computation
// Overloads runSimulation + additional protected methods
#pragma once

#include "computation.h"
#include "pressure_solver/pressureSolverParallel.h"
#include "pressure_solver/sorParallel.h"
#include "output_writer/output_writer_paraview_parallel.h"
#include "output_writer/output_writer_text_parallel.h"

class ComputationParallel : public Computation
{
public:
    /// @brief solves the Navier-Stokes equations by partitioning the domain and therefore calculating parallelly
    ///        the algorithm loops over time and calculates for each time step the updated velocities u and v, the
    ///        pressure p based on the checker board scheme, the preliminary velocities F and G and the right hand side of the pressure Poisson problem (rhs)
    void runSimulation();
   
    /// @brief initializes all objects needed in order to solve the Navier-Stokes equations parallelly
    ///        reads in the parameters from the given parameter file
    ///        decides based on parameter file which discretization scheme to choose
    /// @param argc number of arguments parsed in the command line, must be equal to 2
    /// @param argv the arguments parsed in the command line, contains the file name of the parameter file
    void initialize(int argc, char *argv[]);

private:
    /// @brief if current rank contains Dirichlet boundaries, boundary conditions of u and v
    ///        are set on the local boundary based on given Dirichlet conditions
    ///        handles priority of boundary conditions by introducing shifting parameters for v
    void applyBCOnDirichletBoundary();

    /// @brief if current rank contains Dirichlet boundary, boundary conditions of F and G
    ///        are set on the local boundary based on given Dirichlet conditions
    void applyPreliminaryBCOnDirichletBoundary();

    /// @brief if current rank contains Dirichlet boundary, boundary conditions of u and v
    ///        are set in halo cells based on given Dirichlet conditions and inner cell values
    ///        handles priority of boundary conditions by introducing shifting parameters for u
    void applyBCInHaloCellsAtDirichletBoundary();

    /// @brief diagonal pressure communication for interpolate function and output 
    ///        if the current rank has a lower-left neighbour, the own p(1,1)-value is sent to that neighbour
    ///        if the current rank has an upper-right neighbour, the current rank receives the p(1,1)-value from that neighbour
    void receiveAndSendDiagonalPressureFromAndToOtherProcess();

    /// @brief u and v communication
    ///        checks if current rank contains Dirichlet boundary
    ///        if current rank does not contain a certain Dirichlet boundary, then the calculated u and v need to be communicated with the corresponding neighbour
    ///        function also computes diagonal communication of u and v needed for donor cell scheme
    void receiveAndSendVelocitiesFromAndToOtherProcesses();

    /// @brief computes common time step width for all processes by using serial computeTimeStepWidth()-function
    ///        finds the minimum of all possible time step widths via MPI-function MPI_Allreduce and the operation MPI_MIN
    void computeTimeStepWidthParallel();

    /// @brief computes preliminary velocities F and G
    ///        if current rank does not have a right Dirichlet boundary, then F is also calculated on the right boundary of that rank
    ///        if current rank does not have a top Dirichlet boundary, then G is also calculated on the upper boundary of that rank
    ///        otherwise it has been already set on the Dirichlet boundary by applyPreliminaryBCOnDirichletBoundary() once in the beginning since they don't change
    void computePreliminaryVelocities();

    /// @brief F and G communication
    ///        if current rank does not contain top and/or right Dirichlet boundaries, then G and/or F need to be sent to those neighbours
    ///        if current rank does not contain bottom and/or left Dirichlet boundaries, then G and/or F need to be received from those neighbours
    void receiveAndSendPreliminaryVelocitiesFromAndToOtherProcesses();

    /// @brief computes pressure by using the selected parallel pressure solver SOR-algorithm with checker board pattern
    void computePressure();

    /// @brief computes new velocities u and v
    ///        depending on whether current rank contains Dirichlet boundary or not, the new velocities u and v are also computed on the local boundary of that rank
    ///        this is handled by parameters called shiftIEndU and shiftJEndV
    void computeVelocities();
    
    /// @brief shared pointer to the parallel pressure solver (SOR)
    std::unique_ptr<PressureSolverParallel> pressureSolverParallel_;

    /// @brief unique pointer to the parallel output writer which produces vtk-output
    std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaviewParallel_;

    /// @brief unique pointer to the parallel output writer which produces readable txt-files
    std::unique_ptr<OutputWriterTextParallel> outputWriterTextParallel_;

    /// @brief contains information about the partition of the physical domain needed for MPI
    Partitioning partitioning_;

    /// @brief number of elements in x direction for this process (halo cells not included)
    int nCellsX_;

    /// @brief number of elements in y direction for this process (halo cells not included)
    int nCellsY_;
};