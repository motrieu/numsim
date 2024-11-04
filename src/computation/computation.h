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

    void initialize(int argc, char *argv[]);
 
    void runSimulation();

private:

    void applyBoundaryValues();

    void applyPreliminaryBoundaryValues();
    
    void computeTimeStepWidth();
 	
    void computePreliminaryVelocities();
 	
    void computeRightHandSide();

    void computePressure();
 	
    void computeVelocities();

    Settings settings_;
 
    std::shared_ptr<Discretization> discretization_;
    
    std::unique_ptr<PressureSolver> pressureSolver_;
    
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    
    std::unique_ptr<OutputWriterText> outputWriterText_;
    
    std::array<double,2> meshWidth_;
    
    double dt_;

};
