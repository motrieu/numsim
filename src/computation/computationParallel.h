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
    void runSimulation();
    void initialize(int argc, char *argv[]);

private:
    void applyBCOnDirichletBoundary();
    void applyPreliminaryBCOnDirichletBoundary();
    void applyBCInHaloCellsAtDirichletBoundary();

    //for output
    void receiveAndSendDiagonalPressureFromAndToOtherProcess();
    void receiveAndSendVelocitiesFromAndToOtherProcesses();
    void computeTimeStepWidthParallel();
    void computePreliminaryVelocities();
    void receiveAndSendPreliminaryVelocitiesFromAndToOtherProcesses();
    void computePressure();
    void computeVelocities();
    

    std::unique_ptr<PressureSolverParallel> pressureSolverParallel_;

    std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaviewParallel_;

    std::unique_ptr<OutputWriterTextParallel> outputWriterTextParallel_;

    Partitioning partitioning_;

    /// @brief number of elements in x direction for this process (halo cells not included)
    int nCellsX_;

    /// @brief number of elements in y direction for this process (halo cells not included)
    int nCellsY_;
};