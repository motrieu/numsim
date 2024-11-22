#include "computationParallel.h"


void ComputationParallel::runSimulation();
void ComputationParallel::initialize()
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

    partitioning_ = Partitioning();
    partitioning_.initialize(settings.nCells);

    // either the Central Differences or the Donor cell scheme is used
    if (settings_.useDonorCell)
        discretization_ = std::make_shared<DonorCell>(partitioning_.nCellsLocal(), meshWidth_, settings_.alpha);
    else
        discretization_ = std::make_shared<CentralDifferences>(partitioning_.nCellsLocal(), meshWidth_);

    // either the Gauss-Seidel or the SOR algorithm is used
    if (settings_.pressureSolver == "SOR")
        pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    else if (settings_.pressureSolver == "GaussSeidel")
        pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    else
        throw std::invalid_argument("Only SOR and GaussSeidel are supported as pressure solvers.");
    
    applyBCOnDirichletBoundary();
    applyPreliminaryBCOnDirichletBoundary();
}


void ComputationParallel::applyBCOnDirichletBoundary()
{
    if (partitioning_.ownPartitionContainsLeftBoundary())
    {
        // sets boundary conditions for u(0_local,j_local) based on given Dirichlet conditions
        for (int j=(*discretization_).uJBegin()-1; j < (*discretization_).uJEnd()+1; j++)
            (*discretization_).u((*discretization_).uIBegin()-1,j) = settings_.dirichletBcLeft[0];
    }

    if (partitioning_.ownPartitionContainsRightBoundary())
    {
        // sets boundary conditions for u(N_local,j_local) based on given Dirichlet conditions
        for (int j=(*discretization_).uJBegin()-1; j < (*discretization_).uJEnd()+1; j++)
            (*discretization_).u((*discretization_).uIEnd(),j) = settings_.dirichletBcRight[0];
    }

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        // sets boundary conditions for v(i_local,0_local) based on given Dirichlet conditions
        for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
            (*discretization_).v(i,(*discretization_).vJBegin()-1) = settings_.dirichletBcBottom[1];
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        // sets boundary conditions for v(i_local,N_local) based on given Dirichlet conditions
        for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
            (*discretization_).v(i,(*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
    }
}


void ComputationParallel::applyPreliminaryBCOnDirichletBoundary()
{
    if (partitioning_.ownPartitionContainsLeftBoundary())
    {
        // sets boundary conditions for F(0_local,j_local) based on given Dirichlet conditions
        for (int j=(*discretization_).uJBegin()-1; j < (*discretization_).uJEnd()+1; j++)
            (*discretization_).f((*discretization_).uIBegin()-1,j) = settings_.dirichletBcLeft[0];
    }

    if (partitioning_.ownPartitionContainsRightBoundary())
    {
        // sets boundary conditions for F(N_local,j_local) based on given Dirichlet conditions
        for (int j=(*discretization_).uJBegin()-1; j < (*discretization_).uJEnd()+1; j++)
            (*discretization_).f((*discretization_).uIEnd(),j) = settings_.dirichletBcRight[0];
    }

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        // sets boundary conditions for G(i_local,0_local) based on given Dirichlet conditions
        for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
            (*discretization_).g(i,(*discretization_).vJBegin()-1) = settings_.dirichletBcBottom[1];
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        // sets boundary conditions for G(i_local,N_local) based on given Dirichlet conditions
        for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
            (*discretization_).g(i,(*discretization_).vJEnd()) = settings_.dirichletBcTop[1];
    }
}


void ComputationParallel::applyBCInHaloCellsAtDirichletBoundary()
{
    if (partitioning_.ownPartitionContainsLeftBoundary())
    {
        // sets boundary conditions for v(0_local,j_local) based on given Dirichlet conditions and the inner cell values v(1,j), v(N,j)
        for (int j=(*discretization_).vJBegin()-1; j < (*discretization_).vJEnd()+1; j++)
        {
            const double vLeft = (*discretization_).v((*discretization_).vIBegin(),j);
            (*discretization_).v((*discretization_).vIBegin()-1,j) = 2.0*settings_.dirichletBcLeft[1] - vLeft;  
        }
    }

    if (partitioning_.ownPartitionContainsRightBoundary())
    {
        // sets boundary conditions for v(N+1_local,j_local) based on given Dirichlet conditions and the inner cell values v(1,j), v(N,j)
        for (int j=(*discretization_).vJBegin()-1; j < (*discretization_).vJEnd()+1; j++)
        {
            const double vRight = (*discretization_).v((*discretization_).vIEnd()-1,j);
            (*discretization_).v((*discretization_).vIEnd(),j) = 2.0*settings_.dirichletBcRight[1] - vRight;
        }
    }

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        // sets boundary conditions for u(i_local,0_local) based on given Dirichlet conditions and the inner cell values u(i,1), u(i,N)
        for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
        {
            const double uLower = (*discretization_).u(i,(*discretization_).uJBegin());
            (*discretization_).u(i,(*discretization_).uJBegin()-1) = 2.0*settings_.dirichletBcBottom[0] - uLower;
        }  
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        // sets boundary conditions for u(i_local,N+1_local) based on given Dirichlet conditions and the inner cell values u(i,1), u(i,N)
        for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd(); i++)
        {
            const double uUpper = (*discretization_).u(i,(*discretization_).uJEnd()-1);
            (*discretization_).u(i,(*discretization_).uJEnd()) = 2.0*settings_.dirichletBcTop[0] - uUpper;
        }
    }    
}

void ComputationParallel::applyBCOnNonDirichletBoundary();

void ComputationParallel::applyBCInHaloCellsAtNonDirichletBoundary();
void ComputationParallel::applyPreliminaryBCOnNonDirichletBoundary();
void ComputationParallel::computeTimeStepWidth();