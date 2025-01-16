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

    // calculates mesh width in x- and y-direction based on given parameters
    const double meshWidthX = settings_.physicalSize[0] / settings_.nCells[0];
    const double meshWidthY = settings_.physicalSize[1] / settings_.nCells[1];
    meshWidth_ = {meshWidthX, meshWidthY};

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

    
    bc_ = BoundaryConditions();
    
    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

    // boundary conditions for u and v on the boundary faces only need to be set once in the beginning of the computation
    applyBCOnBoundary();

    // boundary conditions for F and G on the boundary faces only need to be set once in the beginning of the computation
    applyPreliminaryBCOnBoundary();
}

void Computation::applyBoundaryConditions()
{
    for (int i = bc_.setupIBegin(); i < bc_.setupIEnd(); i++)
    {
        for (int j = bc_.setupJBegin(); j < bc_.setupJEnd(); j++)
        {
            if (bc_.setup(i,j) != bc_.indexFluid())
            {
                std::vector<int> edgeDirections;
                if ((j+1 < bc_.setupJEnd()) && (bc_.setup(i,j+1) == bc_.indexFluid()))
                    edgeDirections.push_back(0);
                if ((i+1 < bc_.setupIEnd()) && (bc_.setup(i+1,j) == bc_.indexFluid()))
                    edgeDirections.push_back(1);
                if ((j-1 >= bc_.setupJBegin()) && (bc_.setup(i,j-1) == bc_.indexFluid()))
                    edgeDirections.push_back(2);
                if ((i-1 >= bc_.setupIBegin()) && (bc_.setup(i-1,j) == bc_.indexFluid()))
                    edgeDirections.push_back(3);
                
                if (edgeDirections.size() > 2)
                    throw std::invalid_state("Only corners or edges allowed for obstacles, more than 2 edges given.");
                
                if (edgeDirections.size() == 2)
                {
                    if (edgeDirections[1] != (edgeDirections[0]+1)%4)
                        throw std::invalid_state("Only corners or edges allowed for obstacles, 2 opposite edges given.");

                    if (bc_.setup(i,j) == bc_.indexNoSlip())
                    {
                        bc_.noSlipCorner(i, j, edgeDirections[0]);
                        bc_.pressureNeumannZeroCorner(i, j, edgeDirections[0]);
                    }
                    else
                        throw std::invalid_state("Only NOSLIP allowed for corners of obstacles.");
                }

                else if (edgeDirections.size() == 1)
                {
                    if (bc_.setup(i,j) == bc_.indexNoSlip())
                    {
                        bc_.noSlip(i, j, edgeDirections[0]);
                        bc_.pressureNeumannZero(i, j, edgeDirections[0]);
                    }
                    else if (bc_.setup(i,j) == bc_.indexSlip())
                    {
                        bc_.slip(i, j, edgeDirections[0]);
                        bc_.pressureNeumannZero(i, j, edgeDirections[0]);
                    }
                    else if (bc_.setup(i,j) == bc_.indexInflow())
                    {
                        bc_.inflow(i, j, edgeDirections[0]);
                        bc_.pressureNeumannZero(i, j, edgeDirections[0]);
                    }
                    else if (bc_.setup(i,j) == bc_.indexOutflow())
                    {
                        bc_.outflow(i, j, edgeDirections[0]);
                        bc_.pressureDirichlet(i, j, edgeDirections[0]);
                    }
                }
            }
        }
    }

}

void Computation::computeTimeStepWidth()
{
    const double dx = meshWidth_[0];
    const double dy = meshWidth_[1];
    const double dtDiffusive = (settings_.re/2.0) * (dx*dx * dy*dy) / (dx*dx + dy*dy);
    
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

    const double dtConvectiveU = dx / uAbsMax;
    const double dtConvectiveV = dy / vAbsMax;

    // makes sure that all stability conditions (the convective conditions and the diffusive condition) are fulfilled
    // and that the demanded maximal time step is not exceeded
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