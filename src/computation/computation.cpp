#include "computation.h"

#include <cassert>
#include <array>
#include <memory>
#include <cmath>


void Computation::runSimulation()
{
    double time = 0.0;

    while (time < settings_.endTime)
    {
        // boundary conditions of u and v in halo cells need to be set in each time step
        applyBCInHaloCells();

        // time step width needs to be calculated each time step to ensure stability
        computeTimeStepWidth();

        // ensures that the last time step leads exactly to the demanded end time
        if (time+dt_ > settings_.endTime - dt_/100000.0)
            dt_ = settings_.endTime - time;

        computePreliminaryVelocities();
        
        computeRightHandSide();

        computePressure();

        computeVelocities();

        time += dt_;

        // (*outputWriterParaview_).writeFile(time);
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

    // calculates mesh width in x- and y-direction based on given parameters
    dx_ = settings_.physicalSize[0] / settings_.nCells[0];
    dy_ = settings_.physicalSize[1] / settings_.nCells[1];
    meshWidth_ = {dx_, dy_};
    dxSquared_ = dx_ * dx_;
    dySquared_ = dy_ * dy_;

    // either the Central Differences or the Donor cell scheme is used
    if (settings_.useDonorCell)
        discretization_ = std::make_shared<DonorCell>(settings_.nCells, meshWidth_, settings_.alpha);
    else
        discretization_ = std::make_shared<CentralDifferences>(settings_.nCells, meshWidth_);

    // either the Gauss-Seidel or the SOR algorithm is used
    if (settings_.pressureSolver == "SOR")
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    else if (settings_.pressureSolver == "GaussSeidel")
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    else
        throw std::invalid_argument("Only SOR and GaussSeidel are supported as pressure solvers.");
    
    // outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

    // boundary conditions for u and v on the boundary faces only need to be set once in the beginning of the computation
    applyBCOnBoundary();

    // boundary conditions for F and G on the boundary faces only need to be set once in the beginning of the computation
    applyPreliminaryBCOnBoundary();
}

void Computation::applyBCOnBoundary()
{
    // sets boundary conditions for u(0,j) and u(N,j) based on given Dirichlet conditions
    for (int j=(*discretization_).uJBegin()-1; j < (*discretization_).uJEnd()+1; j++)
    {
        (*discretization_).u((*discretization_).uIBegin()-1,j) = settings_.dirichletBcLeft[0];
        (*discretization_).u((*discretization_).uIEnd(),j) = settings_.dirichletBcRight[0];
    }

    // sets boundary conditions for v(i,0) and v(i,N) based on given Dirichlet conditions
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        (*discretization_).v(i,(*discretization_).vJBegin()-1) = settings_.dirichletBcBottom[1];
        (*discretization_).v(i,(*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
    }
}

void Computation::applyBCInHaloCells()
{
    // sets boundary conditions for u(i,0) and u(i,N+1) based on given Dirichlet conditions and the inner cell values u(i,1), u(i,N)
    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
    {
        const double uLower = (*discretization_).u(i,(*discretization_).uJBegin());
        const double uUpper = (*discretization_).u(i,(*discretization_).uJEnd()-1);
        (*discretization_).u(i,(*discretization_).uJBegin()-1) = 2.0*settings_.dirichletBcBottom[0] - uLower;
        (*discretization_).u(i,(*discretization_).uJEnd()) = 2.0*settings_.dirichletBcTop[0] - uUpper;
    }

    // sets boundary conditions for v(0,j) and v(N+1,j) based on given Dirichlet conditions and the inner cell values v(1,j), v(N,j)
    for (int j=(*discretization_).vJBegin()-1; j < (*discretization_).vJEnd()+1; j++)
    {
        const double vLeft = (*discretization_).v((*discretization_).vIBegin(),j);
        const double vRight = (*discretization_).v((*discretization_).vIEnd()-1,j);
        (*discretization_).v((*discretization_).vIBegin()-1,j) = 2.0*settings_.dirichletBcLeft[1] - vLeft;
        (*discretization_).v((*discretization_).vIEnd(),j) = 2.0*settings_.dirichletBcRight[1] - vRight;
    }
}

void Computation::applyPreliminaryBCOnBoundary()
{
    // sets boundary values of F equal to boundary values of u(0,j), u(N,j)
    for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
    {
        const double uLeft = (*discretization_).u((*discretization_).uIBegin()-1,j);
        const double uRight = (*discretization_).u((*discretization_).uIEnd(),j);
        (*discretization_).f((*discretization_).uIBegin()-1,j) = uLeft;
        (*discretization_).f((*discretization_).uIEnd(),j) = uRight;
    }

    // sets boundary values of G equal to boundary values of v(i,0), v(i,N)
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        const double vLower = (*discretization_).v(i,(*discretization_).vJBegin()-1);
        const double vUpper = (*discretization_).v(i,(*discretization_).vJEnd());
        (*discretization_).g(i,(*discretization_).vJBegin()-1) = vLower;
        (*discretization_).g(i,(*discretization_).vJEnd()) = vUpper;
    }
}

void Computation::computeTimeStepWidth()
{
    const double dtDiffusive = (settings_.re/2.0) * (dxSquared_ * dySquared_) / (dxSquared_ + dySquared_);
    
    double uAbsMax = 0.0;
    for (int i=(*discretization_).uIBegin()-1; i < (*discretization_).uIEnd()+1; i++)
    {
        for (int j=(*discretization_).uJBegin()-1; j < (*discretization_).uJEnd()+1; j++)
        {
            const double uAbs = std::fabs((*discretization_).u(i,j));
            if (uAbs > uAbsMax)
                uAbsMax = uAbs;
        }
    }
    double vAbsMax = 0.0;
    for (int i=(*discretization_).vIBegin()-1; i < (*discretization_).vIEnd()+1; i++)
    {
        for (int j=(*discretization_).vJBegin()-1; j < (*discretization_).vJEnd()+1; j++)
        {
            const double vAbs = std::fabs((*discretization_).v(i,j));
            if (vAbs > vAbsMax)
                vAbsMax = vAbs;
        }
    }

    const double dtConvectiveU = dx_ / uAbsMax;
    const double dtConvectiveV = dy_ / vAbsMax;

    // makes sure that all stability conditions (the convective conditions and the diffusive condition) are fulfilled
    // and that the demanded maximal time step is not exceeded
    dt_ = settings_.tau * std::min({dtDiffusive, dtConvectiveU, dtConvectiveV});
    dt_ = std::min(dt_, settings_.maximumDt);
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
            const double fDiffQuotient = ((*discretization_).f(i,j) - (*discretization_).f(i-1,j)) / dx_;
            const double gDiffQuotient = ((*discretization_).g(i,j) - (*discretization_).g(i,j-1)) / dy_;
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
            (*discretization_).u(i,j) = (*discretization_).f(i,j) - (dt_/dx_)
                                            * ((*discretization_).p(i+1,j) - (*discretization_).p(i,j));
        }
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
        {
            (*discretization_).v(i,j) = (*discretization_).g(i,j) - (dt_/dy_)
                                            * ((*discretization_).p(i,j+1) - (*discretization_).p(i,j));
        }
    }
}