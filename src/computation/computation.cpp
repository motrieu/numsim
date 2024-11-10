#include "computation.h"

#include <cassert>
#include <array>
#include <memory>
#include <cmath>


void Computation::runSimulation()
{
    double time = 0.0;

    // (*outputWriterParaview_).writeFile(time);

    while (time < settings_.endTime)
    {
        computeTimeStepWidth();

        if (time+dt_ > settings_.endTime - dt_/100.0)
            dt_ = settings_.endTime - time;

        computePreliminaryVelocities();
        
        computeRightHandSide();

        computePressure();

        computeVelocities();

        applyBCInHaloCells();

        time += dt_;

        (*outputWriterParaview_).writeFile(time);
        // (*outputWriterText_).writeFile(time);
    }
}

void Computation::initialize(int argc, char *argv[])
{
    assert(argc == 2);

    // read in the first argument
    std::string filename = argv[1];

    // load settings from file
    settings_.loadFromFile(filename);

    const double meshWidthX = settings_.physicalSize[0]/(settings_.nCells[0]-2);
    const double meshWidthY = settings_.physicalSize[1]/(settings_.nCells[1]-2);
    meshWidth_ = {meshWidthX, meshWidthY};

    if (settings_.useDonorCell)
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    else
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);

    if (settings_.pressureSolver == "SOR")
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    else if (settings_.pressureSolver == "GaussSeidel")
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    else
        throw std::invalid_argument("Only SOR and GaussSeidel are supported as pressure solvers.");
    
    applyBCOnBoundary();
    applyBCInHaloCells();
    applyPreliminaryBCOnBoundary();
    // applyPreliminaryBCInHaloCells();

    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);
}

void Computation::applyBCOnBoundary()
{
    for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
    {
        (*discretization_).u((*discretization_).uIBegin()-1,j) = settings_.dirichletBcLeft[0];
        (*discretization_).u((*discretization_).uIEnd(),j) = settings_.dirichletBcRight[0];
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        (*discretization_).v(i,(*discretization_).vJBegin()-1) = settings_.dirichletBcBottom[1];
        (*discretization_).v(i,(*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
    }
}

void Computation::applyBCInHaloCells()
{
    for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
    {
        const double vLeft = (*discretization_).v((*discretization_).vIBegin(),j);
        const double vRight = (*discretization_).v((*discretization_).vIEnd()-1,j);
        (*discretization_).v((*discretization_).vIBegin()-1,j) = 2.0*settings_.dirichletBcLeft[1] - vLeft;
        (*discretization_).v((*discretization_).vIEnd(),j) = 2.0*settings_.dirichletBcRight[1] - vRight;
    }
    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
    {
        const double uLower = (*discretization_).u(i,(*discretization_).uJBegin());
        const double uUpper = (*discretization_).u(i,(*discretization_).uJEnd()-1);
        (*discretization_).u(i,(*discretization_).uJBegin()-1) = 2.0*settings_.dirichletBcBottom[0] - uLower;
        (*discretization_).u(i,(*discretization_).uJEnd()) = 2.0*settings_.dirichletBcTop[0] - uUpper;
    }

    // applyPreliminaryBCInHaloCells();
}

void Computation::applyPreliminaryBCOnBoundary()
{
    for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
    {
        const double uLeft = (*discretization_).u((*discretization_).uIBegin()-1,j);
        const double uRight = (*discretization_).u((*discretization_).uIEnd(),j);
        (*discretization_).f((*discretization_).uIBegin()-1,j) = uLeft;
        (*discretization_).f((*discretization_).uIEnd(),j) = uRight;
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        const double vLower = (*discretization_).v(i,(*discretization_).vJBegin()-1);
        const double vUpper = (*discretization_).v(i,(*discretization_).vJEnd());
        (*discretization_).g(i,(*discretization_).vJBegin()-1) = vLower;
        (*discretization_).g(i,(*discretization_).vJEnd()) = vUpper;
    }
}

/*void Computation::applyPreliminaryBCInHaloCells()
{
    for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
    {
        const double vLeft = (*discretization_).v((*discretization_).vIBegin()-1,j);
        const double vRight = (*discretization_).v((*discretization_).vIEnd(),j);
        (*discretization_).g((*discretization_).vIBegin()-1,j) = vLeft;
        (*discretization_).g((*discretization_).vIEnd(),j) = vRight;
    }
    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
    {
        const double uLower = (*discretization_).u(i,(*discretization_).uJBegin()-1);
        const double uUpper = (*discretization_).u(i,(*discretization_).uJEnd());
        (*discretization_).f(i,(*discretization_).uJBegin()-1) = uLower;
        (*discretization_).f(i,(*discretization_).uJEnd()) = uUpper;
    }
}*/

void Computation::computeTimeStepWidth()
{
    const double dtDiffusive = (settings_.re/2.0) * (meshWidth_[0]*meshWidth_[0]*meshWidth_[1]*meshWidth_[1])
                                        / (meshWidth_[0]*meshWidth_[0] + meshWidth_[1]*meshWidth_[1]);
    
    double uAbsMax = 0.0;
    double vAbsMax = 0.0;
    for (int i=0; i < (*discretization_).nCells()[0]; i++)
    {
        for (int j=0; j < (*discretization_).nCells()[1]; j++)
        {
            const double uAbs = std::fabs((*discretization_).u(i,j));
            const double vAbs = std::fabs((*discretization_).v(i,j));
            if (uAbs > uAbsMax)
                uAbsMax = uAbs;
            if (vAbs > vAbsMax)
                vAbsMax = vAbs;
        }
    }

    const double dtConvectiveU = meshWidth_[0] / uAbsMax;
    const double dtConvectiveV = meshWidth_[1] / vAbsMax;

    dt_ = settings_.tau * std::min({dtDiffusive, dtConvectiveU, dtConvectiveV, settings_.maximumDt});
}

void Computation::computePreliminaryVelocities()
{
    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
    {
        for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
        {
            (*discretization_).f(i,j) = (*discretization_).u(i,j) + dt_ * (
                                            (1.0/settings_.re) * ((*discretization_).computeD2uDx2(i,j) + (*discretization_).computeD2uDy2(i,j))
                                            - (*discretization_).computeDu2Dx(i,j)
                                            - (*discretization_).computeDuvDy(i,j)
                                            + settings_.g[0]
                                            );
        }
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
        {
            (*discretization_).g(i,j) = (*discretization_).v(i,j) + dt_ * (
                                            (1.0/settings_.re) * ((*discretization_).computeD2vDx2(i,j) + (*discretization_).computeD2vDy2(i,j))
                                            - (*discretization_).computeDuvDx(i,j)
                                            - (*discretization_).computeDv2Dy(i,j)
                                            + settings_.g[1]
                                            );
        }
    }
}

void Computation::computeRightHandSide()
{
    for (int i=(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
    {
        for (int j=(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
        {
            const double fDiffQuotient = ((*discretization_).f(i,j) - (*discretization_).f(i-1,j)) / meshWidth_[0];
            const double gDiffQuotient = ((*discretization_).g(i,j) - (*discretization_).g(i,j-1)) / meshWidth_[1];
            (*discretization_).rhs(i,j) = (1.0/dt_) * (fDiffQuotient + gDiffQuotient);
        }
    }
}

void Computation::computePressure()
{
    (*pressureSolver_).solve();
}

void Computation::computeVelocities()
{
    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
    {
        for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
        {
            (*discretization_).u(i,j) = (*discretization_).f(i,j) - (dt_/meshWidth_[0])
                                            * ((*discretization_).p(i+1,j) - (*discretization_).p(i,j));
        }
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
        {
            (*discretization_).v(i,j) = (*discretization_).g(i,j) - (dt_/meshWidth_[1])
                                            * ((*discretization_).p(i,j+1) - (*discretization_).p(i,j));
        }
    }
}