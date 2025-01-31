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

    applyBoundaryConditions();
    (*outputWriterParaview_).writeFile(time);
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
    // outputWriterText_ = std::make_unique<OutputWriterText>(discretization_);

    // read in setup from parameter file
    loadSetupFromFile(filename);

    // needed for applying boundary conditions in the right directions
    initializeEdgeDirections();
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
            //edgeDirections are only written in obstacle cells since boundary conditions are applied from obstacle pov
            if ((*discretization_).setup(i,j) != (*discretization_).indexFluid())
            {
                std::vector<int> edgeDirs;
                //if upper neighbour is fluid, index 0 is pushed as north direction
                if ((j+1 < (*discretization_).setupJEnd()) && ((*discretization_).setup(i,j+1) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(0);
                //if right neighbour is fluid, index 1 is pushed as east direction
                if ((i+1 < (*discretization_).setupIEnd()) && ((*discretization_).setup(i+1,j) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(1);
                //if lower neighbour is fluid, index 2 is pushed as south direction
                if ((j-1 >= (*discretization_).setupJBegin()) && ((*discretization_).setup(i,j-1) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(2);
                //if left neighbour is fluid, index 3 is pushed as west direction
                if ((i-1 >= (*discretization_).setupIBegin()) && ((*discretization_).setup(i-1,j) == (*discretization_).indexFluid()))
                    edgeDirs.push_back(3);
                
                //an obstacle cannot have more than two edgeDirections to prevent information from passing through
                if (edgeDirs.size() > 2)
                    throw std::invalid_argument("Only corners or edges allowed for obstacles, more than 2 edges given.");
                
                //if an obstacle has two edges, only a corner is allowed
                if (edgeDirs.size() == 2)
                {
                    int edgeDirsDiff = edgeDirs[1]-edgeDirs[0];

                    //no opposite edges as edgeDirections are allowed to prevent information from passing through
                    if (edgeDirsDiff == 2)
                        throw std::invalid_argument("Only corners or edges allowed for obstacles, 2 opposite edges given.");
                    if ((*discretization_).setup(i,j) == (*discretization_).indexNoSlip())
                    {
                        //if the difference of the indices listed in edgeDirs is 1, the possible corners are 1 for north-east, 3 for east-south and 5 for south-west
                        if (edgeDirsDiff == 1)
                            (*discretization_).edgeDirections(i,j) = edgeDirs[0]*2 + 1;

                        //if the difference of the indices listed in edgeDirs is 3, only the west-north corner is possible and therefore a 7 is stored in edgeDirections
                        else if (edgeDirsDiff == 3)
                            (*discretization_).edgeDirections(i,j) = 7;
                        else
                            throw std::invalid_argument("Invalid edge-direction combination given.");
                    }
                    else
                        throw std::invalid_argument("Only NOSLIP allowed for corners of obstacles.");
                }
                //if an obstacle only has one edge either 0 for north, 2 for east, 4 for south and 6 for west is stored in edgeDirections
                else if (edgeDirs.size() == 1)
                {
                    (*discretization_).edgeDirections(i,j) = edgeDirs[0]*2;
                }
                //if an obstacle has no edges at all, one still needs to check for diagonal fluid cells
                else
                {
                    std::vector<int> diagonalFluidCells;
                    //upper right neighbour is fluid cell, index 1 is pushed
                    if ((i+1 < (*discretization_).setupIEnd())
                            && (j+1 < (*discretization_).setupJEnd())
                            && ((*discretization_).setup(i+1,j+1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(1);
                    
                    //lower right neighbour is fluid cell, index 3 is pushed
                    if ((i+1 < (*discretization_).setupIEnd())
                            && (j-1 >= (*discretization_).setupJBegin())
                            && ((*discretization_).setup(i+1,j-1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(3);
                    
                    //lower left neighbour is fluid cell, index 5 is pushed
                    if ((i-1 >= (*discretization_).setupIBegin())
                            && (j-1 >= (*discretization_).setupJBegin())
                            && ((*discretization_).setup(i-1,j-1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(5);
                    
                    //upper right neighbour is fluid cell, index 7 is pushed
                    if ((i-1 >= (*discretization_).setupIBegin())
                            && (j+1 < (*discretization_).setupJEnd())
                            && ((*discretization_).setup(i-1,j+1) == (*discretization_).indexFluid()))
                        diagonalFluidCells.push_back(7);
                    
                    if (diagonalFluidCells.size() > 1)
                        throw std::invalid_argument("Only one diagonal adjoining fluid cell allowed, more than 1 given.");
                    else if (diagonalFluidCells.size() == 1)
                        (*discretization_).edgeDirections(i,j) = diagonalFluidCells[0];
                }

                //numberFaces stores the number of edges an obstacle has towards fluid cells
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
            //boundary values for obstacle cells are set depending on the value written down in the setup array
            //only obstacle cells with minimum 1 fluid edge are handled in this part
            //no obstacle cells with diagonal fluid cells are handled in this part
            if (((*discretization_).setup(i,j) != (*discretization_).indexFluid())
                    && ((*discretization_).edgeDirections(i,j) != -1)
                    && ((*discretization_).numberFaces(i,j) != 0))
            {
                int edgeDirection = (*discretization_).edgeDirections(i,j);

                //if the edgeDirection indicates a corner the NOSLIP-corner condition for u and v is applied
                if (edgeDirection%2 == 1)
                {
                    (*discretization_).noSlipCorner(i, j, edgeDirection);
                }
                //if the edgeDirection indicates a single edge the boundary condition written down in the setup array is applied for u and v
                //can be either NOSLIP, SLIP, INFLOW, OUTFLOW, PRESSURE
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
                //if numberFaces is zero, then the obstacle has a diagonal fluid cell 
                //this is why then the diagonal boundary functions are used to apply the boundary conditions
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
                //setting g and f boundary conditions depending on edgeDirections (corner or edge)
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
    //applying boundary conditions for pressure
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