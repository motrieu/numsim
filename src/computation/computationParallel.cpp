#include "computationParallel.h"


void ComputationParallel::runSimulation()
{
    receiveAndSendVelocitiesFromAndToOtherProcesses();
    computeTimeStepWidthParallel();
}
void ComputationParallel::initialize(int argc, char *argv[])
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
    partitioning_.initialize(settings_.nCells);

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

void ComputationParallel::receiveAndSendVelocitiesFromAndToOtherProcesses()
{
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> receiveLeftRequests;
    std::vector<MPI_Request> receiveRightRequests;
    std::vector<MPI_Request> receiveLowerRequests;
    std::vector<MPI_Request> receiveUpperRequests;

    int nCellsX = (*discretization_).nCells()[0];
    int nCellsY = (*discretization_).nCells()[1];

    std::vector<double> receiveLeftUBuffer(nCellsY);
    std::vector<double> receiveLeftVBuffer(nCellsY);
    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        std::vector<double> sendLeftUBuffer(nCellsY);
        std::vector<double> sendLeftVBuffer(nCellsY);
       
        for (int j = 1; j < nCellsY+1; j++)
        {
            sendLeftUBuffer[j-1] = (*discretization_).u(1,j);
            sendLeftVBuffer[j-1] = (*discretization_).v(1,j);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendLeftUBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendLeftVBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveLeftRequests.emplace_back();
        MPI_Irecv(receiveLeftUBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLeftRequests.back());

        receiveLeftRequests.emplace_back();
        MPI_Irecv(receiveLeftVBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveLeftRequests.back());
    }

    std::vector<double> receiveRightUBuffer(nCellsY);
    std::vector<double> receiveRightVBuffer(nCellsY);
    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        std::vector<double> sendRightUBuffer(nCellsY);
        std::vector<double> sendRightVBuffer(nCellsY);
       
        for (int j = 1; j < nCellsY+1; j++)
        {
            sendRightUBuffer[j-1] = (*discretization_).u(nCellsX,j);
            sendRightVBuffer[j-1] = (*discretization_).v(nCellsX,j);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendRightUBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendRightVBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveRightRequests.emplace_back();
        MPI_Irecv(receiveRightUBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRightRequests.back());

        receiveRightRequests.emplace_back();
        MPI_Irecv(receiveRightVBuffer.data(), nCellsY, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveRightRequests.back());
    }

    std::vector<double> receiveLowerUBuffer(nCellsX);
    std::vector<double> receiveLowerVBuffer(nCellsX);
    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerUBuffer(nCellsX);
        std::vector<double> sendLowerVBuffer(nCellsX);
       
        for (int i = 1; i < nCellsX+1; i++)
        {
            sendLowerUBuffer[i-1] = (*discretization_).u(i,1);
            sendLowerVBuffer[i-1] = (*discretization_).v(i,1);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendLowerUBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendLowerVBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveLowerRequests.emplace_back();
        MPI_Irecv(receiveLowerUBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLowerRequests.back());

        receiveLowerRequests.emplace_back();
        MPI_Irecv(receiveLowerVBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveLowerRequests.back());
    }

    std::vector<double> receiveUpperUBuffer(nCellsX);
    std::vector<double> receiveUpperVBuffer(nCellsX);
    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendUpperUBuffer(nCellsX);
        std::vector<double> sendUpperVBuffer(nCellsX);
       
        for (int i = 1; i < nCellsX+1; i++)
        {
            sendUpperUBuffer[i-1] = (*discretization_).u(i,nCellsY);
            sendUpperVBuffer[i-1] = (*discretization_).v(i,nCellsY);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendUpperUBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendUpperVBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveUpperRequests.emplace_back();
        MPI_Irecv(receiveUpperUBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveUpperRequests.back());

        receiveUpperRequests.emplace_back();
        MPI_Irecv(receiveUpperVBuffer.data(), nCellsX, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveUpperRequests.back());
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Waitall(receiveLeftRequests.size(), receiveLeftRequests.data(), MPI_STATUSES_IGNORE);
        for (int j = 1; j < nCellsY+1; j++)
        {
            (*discretization_).u(0,j) = receiveLeftUBuffer[j-1];
            (*discretization_).v(0,j) = receiveLeftVBuffer[j-1];
        }  
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Waitall(receiveRightRequests.size(), receiveRightRequests.data(), MPI_STATUSES_IGNORE);
        for (int j = 1; j < nCellsY+1; j++)
        {
            (*discretization_).u(nCellsX+1,j) = receiveRightUBuffer[j-1];
            (*discretization_).v(nCellsX+1,j) = receiveRightVBuffer[j-1];
        }  
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Waitall(receiveLowerRequests.size(), receiveLowerRequests.data(), MPI_STATUSES_IGNORE);
        for (int i = 1; i < nCellsX+1; i++)
        {
            (*discretization_).u(i,0) = receiveLowerUBuffer[i-1];
            (*discretization_).v(i,0) = receiveLowerVBuffer[i-1];
        }  
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Waitall(receiveUpperRequests.size(), receiveUpperRequests.data(), MPI_STATUSES_IGNORE);
        for (int i = 1; i < nCellsX+1; i++)
        {
            (*discretization_).u(i,nCellsY+1) = receiveUpperUBuffer[i-1];
            (*discretization_).v(i,nCellsY+1) = receiveUpperVBuffer[i-1];
        }  
    }

    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
}

void ComputationParallel::computeTimeStepWidthParallel()
{
    computeTimeStepWidth();

    const double sendDt = dt_;
    double receiveDt;
    
    MPI_Allreduce(&sendDt, &receiveDt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    dt_ = receiveDt;
}

//use offset for pressure solver checkerboard pattern