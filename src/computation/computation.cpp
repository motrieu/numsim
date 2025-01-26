#include "computation.h"

#include <cassert>
#include <array>
#include <memory>
#include <cmath>

#include <iostream>


void Computation::runSimulation()
{
    double time = 0.0;

    while (time < settings_.endTime)
    {
        std::cout << time << std::endl;

        // boundary conditions of u and v in halo cells need to be set in each time step
        applyBoundaryConditions();

        (*outputWriterParaview_).writeFile(time);

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
    }
}

void Computation::initialize(int argc, char *argv[])
{
    assert(argc == 2);

    // read in the first argument
    std::string filename = argv[1];

    // load settings from file
    settings_.loadParamsFromFile(filename);

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

    outputWriterParaview_ = std::make_unique<OutputWriterParaview>(discretization_);

    loadSetupFromFile(filename);

    initializeEdgeDirections();

    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

    // boundary conditions for u and v on the boundary faces only need to be set once in the beginning of the computation
    //applyBCOnBoundary();

    // boundary conditions for F and G on the boundary faces only need to be set once in the beginning of the computation
    //applyPreliminaryBCOnBoundary();
}

void Computation::loadSetupFromFile(std::string filename)
{
    // open file
    std::ifstream file(filename.c_str(), std::ios::in);

    // check if file is open
    if (!file.is_open())
    {
        std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
        return;
    }

    bool readSetup = false;
    bool readUIn = false;
    bool readVIn = false;
    bool readPRB = false;

    int j = 0;

    // loop over lines of file
    for (int lineNo = 0;; lineNo++)
    {
        // read line
        std::string line;
        getline(file, line);

        // at the end of the file break for loop
        if (file.eof())
            break;

        settings_.removeWhitespaceAtBeginning(line);

        if (line.find("#SETUP") != std::string::npos)
        {
            j = 0;
            readSetup = true;
        }
        else if (line.find("#UIN") != std::string::npos)
        {
            j = 0;
            readUIn = true;
        }
        else if (line.find("#VIN") != std::string::npos)
        {
            j = 0;
            readVIn = true;
        }
        else if (line.find("#PRB") != std::string::npos)
        {
            j = 0;
            readPRB = true;
        }
        
        if (readPRB && (std::isdigit(line[0]) || std::isdigit(line[1])))
        {
            int valueStartIndex = 0;
            int valueEndIndex = 0;
            std::string valueString;
            int i = 0;
            while (valueStartIndex < line.length())
            {
                valueEndIndex = line.find_first_of(" \n", valueStartIndex);
                valueString = line.substr(valueStartIndex, valueEndIndex-valueStartIndex);
                double value = std::stod(valueString);
                valueStartIndex = valueEndIndex+1;
                (*discretization_).pRB(i,(*discretization_).setupJEnd()-1-j) = value;
                i++;
            }
            j++;
        }
        else if (readVIn && (std::isdigit(line[0]) || std::isdigit(line[1])))
        {
            int valueStartIndex = 0;
            int valueEndIndex = 0;
            std::string valueString;
            int i = 0;
            while (valueStartIndex < line.length())
            {
                valueEndIndex = line.find_first_of(" \n", valueStartIndex);
                valueString = line.substr(valueStartIndex, valueEndIndex-valueStartIndex);
                double value = std::stod(valueString);
                valueStartIndex = valueEndIndex+1;
                (*discretization_).vIn(i,(*discretization_).setupJEnd()-1-j) = value;
                i++;
            }
            j++;
        }
        else if (readUIn && (std::isdigit(line[0]) || std::isdigit(line[1])))
        {
            int valueStartIndex = 0;
            int valueEndIndex = 0;
            std::string valueString;
            int i = 0;
            while (valueStartIndex < line.length())
            {
                valueEndIndex = line.find_first_of(" \n", valueStartIndex);
                valueString = line.substr(valueStartIndex, valueEndIndex-valueStartIndex);
                double value = std::stod(valueString);
                valueStartIndex = valueEndIndex+1;
                (*discretization_).uIn(i,(*discretization_).setupJEnd()-1-j) = value;
                i++;
            }
            j++;
        }
        else if (readSetup && (std::isdigit(line[0]) || std::isdigit(line[1])))
        {
            int valueStartIndex = 0;
            int valueEndIndex = 0;
            std::string valueString;
            int i = 0;
            while (valueStartIndex < line.length())
            {
                valueEndIndex = line.find_first_of(" \n", valueStartIndex);
                valueString = line.substr(valueStartIndex, valueEndIndex-valueStartIndex);
                int value = std::stod(valueString);
                valueStartIndex = valueEndIndex+1;
                (*discretization_).setup(i,(*discretization_).setupJEnd()-1-j) = value;
                i++;
            }
            j++;
        }
    }
}

void Computation::initializeEdgeDirections()
{
    for (int i = (*discretization_).setupIBegin(); i < (*discretization_).setupIEnd(); i++)
    {
        for (int j = (*discretization_).setupJBegin(); j < (*discretization_).setupJEnd(); j++)
        {
            if ((*discretization_).setup(i,j) != (*discretization_).indexFluid())
            {
                std::vector<int> edgeDirs;
                if ((j+1 < (*discretization_).setupJEnd()) && ((*discretization_).setup(i,j+1) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(0);
                if ((i+1 < (*discretization_).setupIEnd()) && ((*discretization_).setup(i+1,j) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(1);
                if ((j-1 >= (*discretization_).setupJBegin()) && ((*discretization_).setup(i,j-1) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(2);
                if ((i-1 >= (*discretization_).setupIBegin()) && ((*discretization_).setup(i-1,j) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(3);
                
                if (edgeDirs.size() > 2)
                    throw std::invalid_argument("Only corners or edges allowed for obstacles, more than 2 edges given.");
                
                if (edgeDirs.size() == 2)
                {
                    if ((edgeDirs[1]-edgeDirs[0]) == 2)
                        throw std::invalid_argument("Only corners or edges allowed for obstacles, 2 opposite edges given.");

                    if ((*discretization_).setup(i,j) == (*discretization_).indexNoSlip())
                        (*discretization_).edgeDirections(i,j) = edgeDirs[0]*2 + 1;
                    else
                        throw std::invalid_argument("Only NOSLIP allowed for corners of obstacles.");
                }
                else if (edgeDirs.size() == 1)
                {
                    (*discretization_).edgeDirections(i,j) = edgeDirs[0]*2;
                }
                else
                {
                    std::vector<int> diagonalFluidCells;
                    if ((i+1 < (*discretization_).setupIEnd())
                            && (j+1 < (*discretization_).setupJEnd())
                            && ((*discretization_).setup(i+1,j+1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(1);
                    if ((i+1 < (*discretization_).setupIEnd())
                            && (j-1 >= (*discretization_).setupJBegin())
                            && ((*discretization_).setup(i+1,j-1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(3);
                    if ((i-1 >= (*discretization_).setupIBegin())
                            && (j-1 >= (*discretization_).setupJBegin())
                            && ((*discretization_).setup(i-1,j-1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(5);
                    if ((i-1 >= (*discretization_).setupIBegin())
                            && (j+1 < (*discretization_).setupJEnd())
                            && ((*discretization_).setup(i-1,j+1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(7);
                    
                    if (diagonalFluidCells.size() > 1)
                        throw std::invalid_argument("Only one diagonal adjoining fluid cell allowed, more than 1 given.");
                    else if (diagonalFluidCells.size() == 1)
                        (*discretization_).edgeDirections(i,j) = diagonalFluidCells[0];
                }

                (*discretization_).numberFaces(i,j) = edgeDirs.size();
            }
        }
    }
}

void Computation::applyBoundaryConditions()
{
    for (int i = (*discretization_).setupIBegin(); i < (*discretization_).setupIEnd(); i++)
    {
        for (int j = (*discretization_).setupJBegin(); j < (*discretization_).setupJEnd(); j++)
        {
            if (((*discretization_).setup(i,j) != (*discretization_).indexFluid())
                    && ((*discretization_).edgeDirections(i,j) != -1)
                    && ((*discretization_).numberFaces(i,j) != 0))
            {
                int edgeDirection = (*discretization_).edgeDirections(i,j);
                if (edgeDirection%2 == 1)
                {
                    (*discretization_).noSlipCorner(i, j, edgeDirection);
                }
                else
                {
                    if ((*discretization_).setup(i,j) == (*discretization_).indexNoSlip())
                        (*discretization_).noSlip(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexSlip())
                        (*discretization_).slip(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexInflow())
                        (*discretization_).inflow(i, j, edgeDirection, (*discretization_).uIn(i,j), (*discretization_).vIn(i,j));
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexOutflow())
                        (*discretization_).outflow(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexPressure())
                        (*discretization_).outflow(i, j, edgeDirection);
                }
            }
        }
    }

    for (int i = (*discretization_).setupIBegin(); i < (*discretization_).setupIEnd(); i++)
    {
        for (int j = (*discretization_).setupJBegin(); j < (*discretization_).setupJEnd(); j++)
        {
            if (((*discretization_).setup(i,j) != (*discretization_).indexFluid())
                    && ((*discretization_).edgeDirections(i,j) != -1))
            {
                int edgeDirection = (*discretization_).edgeDirections(i,j);
                if ((*discretization_).numberFaces(i,j) == 0)
                {
                    if ((*discretization_).setup(i,j) == (*discretization_).indexNoSlip())
                        (*discretization_).noSlipDiagonal(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexSlip())
                        (*discretization_).slipDiagonal(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexInflow())
                        (*discretization_).inflowDiagonal(i, j, edgeDirection, (*discretization_).uIn(i,j), (*discretization_).vIn(i,j));
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexOutflow())
                        (*discretization_).outflowDiagonal(i, j, edgeDirection);
                    else if ((*discretization_).setup(i,j) == (*discretization_).indexPressure())
                        (*discretization_).outflowDiagonal(i, j, edgeDirection);
                }
                else
                {
                    if (edgeDirection == 0)
                    {
                        (*discretization_).g(i,j) = (*discretization_).v(i,j);
                    }
                    else if (edgeDirection == 1)
                    {
                        (*discretization_).f(i,j) = (*discretization_).u(i,j);
                        (*discretization_).g(i,j) = (*discretization_).v(i,j);
                    }
                    else if (edgeDirection == 2)
                    {
                        (*discretization_).f(i,j) = (*discretization_).u(i,j);
                    }
                    else if (edgeDirection == 3)
                    {
                        (*discretization_).f(i,j) = (*discretization_).u(i,j);
                        (*discretization_).g(i,j-1) = (*discretization_).v(i,j-1);
                    }
                    else if (edgeDirection == 4)
                    {
                        (*discretization_).g(i,j-1) = (*discretization_).v(i,j-1);
                    }
                    else if (edgeDirection == 5)
                    {
                        (*discretization_).f(i-1,j) = (*discretization_).u(i-1,j);
                        (*discretization_).g(i,j-1) = (*discretization_).v(i,j-1);
                    }
                    else if (edgeDirection == 6)
                    {
                        (*discretization_).f(i-1,j) = (*discretization_).u(i-1,j);
                    }
                    else if (edgeDirection == 7)
                    {
                        (*discretization_).f(i-1,j) = (*discretization_).u(i-1,j);
                        (*discretization_).g(i,j) = (*discretization_).v(i,j);
                    }
                }
            }
        }
    }

    (*pressureSolver_).applyBoundaryConditions();
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
    for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
    {
        for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
        {
            if (((*discretization_).setup(i,j) == (*discretization_).indexFluid())
                    && ((*discretization_).setup(i+1,j) == (*discretization_).indexFluid()))
            {
                (*discretization_).f(i,j) = (*discretization_).u(i,j) + dt_ * (
                                            (1.0/settings_.re) * ((*discretization_).computeD2uDx2(i,j) + (*discretization_).computeD2uDy2(i,j))
                                            - (*discretization_).computeDu2Dx(i,j)
                                            - (*discretization_).computeDuvDy(i,j)
                                            + settings_.g[0]
                                            );
            }
        }
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
        {
            if (((*discretization_).setup(i,j) == (*discretization_).indexFluid())
                    && ((*discretization_).setup(i,j+1) == (*discretization_).indexFluid()))
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
}

void Computation::computeRightHandSide()
{
    for (int i=(*discretization_).pIBegin(); i < (*discretization_).pIEnd(); i++)
    {
        for (int j=(*discretization_).pJBegin(); j < (*discretization_).pJEnd(); j++)
        {
            if ((*discretization_).setup(i,j) == (*discretization_).indexFluid())
            {
                const double fDiffQuotient = ((*discretization_).f(i,j) - (*discretization_).f(i-1,j)) / meshWidth_[0];
                const double gDiffQuotient = ((*discretization_).g(i,j) - (*discretization_).g(i,j-1)) / meshWidth_[1];
                (*discretization_).rhs(i,j) = (1.0/dt_) * (fDiffQuotient + gDiffQuotient);
            }
        }
    }
}

void Computation::computePressure()
{
    (*pressureSolver_).solve();
}

void Computation::computeVelocities()
{
    for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
    {
        for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
        {
            if (((*discretization_).setup(i,j) == (*discretization_).indexFluid())
                    && ((*discretization_).setup(i+1,j) == (*discretization_).indexFluid()))
            {
                (*discretization_).u(i,j) = (*discretization_).f(i,j) - (dt_/meshWidth_[0])
                                            * ((*discretization_).p(i+1,j) - (*discretization_).p(i,j));
            }
        }
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd(); j++)
        {
            if (((*discretization_).setup(i,j) == (*discretization_).indexFluid())
                    && ((*discretization_).setup(i,j+1) == (*discretization_).indexFluid()))
            {
                (*discretization_).v(i,j) = (*discretization_).g(i,j) - (dt_/meshWidth_[1])
                                            * ((*discretization_).p(i,j+1) - (*discretization_).p(i,j));
            }
        }
    }
}