#include "computationParallel.h"

#include <mpi.h>

void ComputationParallel::runSimulation()
{
    double time = 0.0;
    const double outputIntervall = 1.0;
    double timeNextOutput = outputIntervall;

    while (time < settings_.endTime)
    {
        receiveAndSendVelocitiesFromAndToOtherProcesses();

        if (time >= timeNextOutput)
        {
            receiveAndSendDiagonalPressureFromAndToOtherProcess();

            (*outputWriterParaviewParallel_).writeFile(time);
            //(*outputWriterTextParallel_).writeFile(time);
            timeNextOutput += outputIntervall;
        }

        applyBCInHaloCellsAtDirichletBoundary();

        computeTimeStepWidthParallel();

        // ensures that the last time step leads exactly to the demanded end time
        if (time+dt_ > settings_.endTime - dt_/100000.0)
            dt_ = settings_.endTime - time;

        computePreliminaryVelocities();

        receiveAndSendPreliminaryVelocitiesFromAndToOtherProcesses();

        computeRightHandSide();

        computePressure();

        computeVelocities();

        time += dt_;
    }
    
    receiveAndSendVelocitiesFromAndToOtherProcesses();
    receiveAndSendDiagonalPressureFromAndToOtherProcess();
    (*outputWriterParaviewParallel_).writeFile(time);
}
void ComputationParallel::initialize(int argc, char *argv[])
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

    partitioning_ = Partitioning();
    partitioning_.initialize(settings_.nCells);

    // either the Central Differences or the Donor cell scheme is used
    if (settings_.useDonorCell)
        discretization_ = std::make_shared<DonorCell>(partitioning_.nCellsLocal(), meshWidth_, settings_.alpha);
    else
        discretization_ = std::make_shared<CentralDifferences>(partitioning_.nCellsLocal(), meshWidth_);

    // either the Gauss-Seidel or the SOR algorithm is used
    if (settings_.pressureSolver == "SOR")
        pressureSolverParallel_ = std::make_unique<SORParallel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_, settings_.omega);
    //else if (settings_.pressureSolver == "GaussSeidel")
    //    pressureSolver_ = std::make_unique<GaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations);
    else
        throw std::invalid_argument("Only SOR and GaussSeidel are supported as pressure solvers.");
    
    outputWriterParaviewParallel_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, partitioning_);
    outputWriterTextParallel_ = std::make_unique<OutputWriterTextParallel>(discretization_, partitioning_);

    nCellsX_ = (*discretization_).nCells()[0];
    nCellsY_ = (*discretization_).nCells()[1];

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
        int shiftIBeginV = 0;
        int shiftIEndV = 0;
        if (!partitioning_.ownPartitionContainsLeftBoundary())
            shiftIBeginV = -1;
        if (!partitioning_.ownPartitionContainsRightBoundary())
            shiftIEndV = 1;

        // sets boundary conditions for v(i_local,0_local) based on given Dirichlet conditions
        for (int i=(*discretization_).vIBegin()+shiftIBeginV; i < (*discretization_).vIEnd()+shiftIEndV; i++)
            (*discretization_).v(i,(*discretization_).vJBegin()-1) = settings_.dirichletBcBottom[1];
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        int shiftIBeginV = 0;
        int shiftIEndV = 0;
        if (!partitioning_.ownPartitionContainsLeftBoundary())
            shiftIBeginV = -1;
        if (!partitioning_.ownPartitionContainsRightBoundary())
            shiftIEndV = 1;

        // sets boundary conditions for v(i_local,N_local) based on given Dirichlet conditions
        for (int i=(*discretization_).vIBegin()+shiftIBeginV; i < (*discretization_).vIEnd()+shiftIEndV; i++)
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
        int shiftIBeginU = 0;
        int shiftIEndU = 0;
        if (!partitioning_.ownPartitionContainsLeftBoundary())
            shiftIBeginU = -1;
        if (!partitioning_.ownPartitionContainsRightBoundary())
            shiftIEndU = 1;

        // sets boundary conditions for u(i_local,0_local) based on given Dirichlet conditions and the inner cell values u(i,1), u(i,N)
        for (int i=(*discretization_).uIBegin()+shiftIBeginU; i < (*discretization_).uIEnd()+shiftIEndU; i++)
        {
            const double uLower = (*discretization_).u(i,(*discretization_).uJBegin());
            (*discretization_).u(i,(*discretization_).uJBegin()-1) = 2.0*settings_.dirichletBcBottom[0] - uLower;
        }  
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        int shiftIBeginU = 0;
        int shiftIEndU = 0;
        if (!partitioning_.ownPartitionContainsLeftBoundary())
            shiftIBeginU = -1;
        if (!partitioning_.ownPartitionContainsRightBoundary())
            shiftIEndU = 1;

        // sets boundary conditions for u(i_local,N+1_local) based on given Dirichlet conditions and the inner cell values u(i,1), u(i,N)
        for (int i=(*discretization_).uIBegin()+shiftIBeginU; i < (*discretization_).uIEnd()+shiftIEndU; i++)
        {
            const double uUpper = (*discretization_).u(i,(*discretization_).uJEnd()-1);
            (*discretization_).u(i,(*discretization_).uJEnd()) = 2.0*settings_.dirichletBcTop[0] - uUpper;
        }
    }    
}

void ComputationParallel::receiveAndSendDiagonalPressureFromAndToOtherProcess()
{
    std::vector<MPI_Request> sendLeftLowerDiagonalRequest;
    MPI_Request receiveRightUpperDiagonalRequest;

    if ((!partitioning_.ownPartitionContainsLeftBoundary()) && (!partitioning_.ownPartitionContainsBottomBoundary()))
    {
        std::vector<double> sendLeftLowerDiagonalPBuffer = {(*discretization_).p(1,1)};

        sendLeftLowerDiagonalRequest.emplace_back();
        MPI_Isend(sendLeftLowerDiagonalPBuffer.data(), 1, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo()-1, 0, MPI_COMM_WORLD, &sendLeftLowerDiagonalRequest.back());
    }

    std::vector<double> receiveRightUpperDiagonalPBuffer(1, 0);
    if ((!partitioning_.ownPartitionContainsRightBoundary()) && (!partitioning_.ownPartitionContainsTopBoundary()))
    {
        MPI_Irecv(receiveRightUpperDiagonalPBuffer.data(), 1, MPI_DOUBLE, partitioning_.topNeighbourRankNo()+1, 0, MPI_COMM_WORLD, &receiveRightUpperDiagonalRequest);
    }

    MPI_Waitall(sendLeftLowerDiagonalRequest.size(), sendLeftLowerDiagonalRequest.data(), MPI_STATUSES_IGNORE);

    if ((!partitioning_.ownPartitionContainsRightBoundary()) && (!partitioning_.ownPartitionContainsTopBoundary()))
    {
        MPI_Wait(&receiveRightUpperDiagonalRequest, MPI_STATUS_IGNORE);

        (*discretization_).p(nCellsX_+1,nCellsY_+1) = receiveRightUpperDiagonalPBuffer.at(0);
    }
}

void ComputationParallel::receiveAndSendVelocitiesFromAndToOtherProcesses()
{
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> receiveLeftRequests;
    std::vector<MPI_Request> receiveRightRequests;
    std::vector<MPI_Request> receiveLowerRequests;
    std::vector<MPI_Request> receiveUpperRequests;

    std::vector<double> receiveLeftUBuffer(nCellsY_, 0);
    std::vector<double> receiveLeftVBuffer(nCellsY_, 0);
    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        std::vector<double> sendLeftUBuffer(nCellsY_, 0);
        std::vector<double> sendLeftVBuffer(nCellsY_, 0);
       
        for (int j = 1; j < nCellsY_+1; j++)
        {
            sendLeftUBuffer.at(j-1) = (*discretization_).u(1,j);
            sendLeftVBuffer.at(j-1) = (*discretization_).v(1,j);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendLeftUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendLeftVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveLeftRequests.emplace_back();
        MPI_Irecv(receiveLeftUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLeftRequests.back());

        receiveLeftRequests.emplace_back();
        MPI_Irecv(receiveLeftVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveLeftRequests.back());
    }

    std::vector<double> receiveRightUBuffer(nCellsY_, 0);
    std::vector<double> receiveRightVBuffer(nCellsY_, 0);
    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        std::vector<double> sendRightUBuffer(nCellsY_, 0);
        std::vector<double> sendRightVBuffer(nCellsY_, 0);
       
        for (int j = 1; j < nCellsY_+1; j++)
        {
            sendRightUBuffer.at(j-1) = (*discretization_).u(nCellsX_,j);
            sendRightVBuffer.at(j-1) = (*discretization_).v(nCellsX_,j);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendRightUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendRightVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveRightRequests.emplace_back();
        MPI_Irecv(receiveRightUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRightRequests.back());

        receiveRightRequests.emplace_back();
        MPI_Irecv(receiveRightVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveRightRequests.back());
    }

    std::vector<double> receiveLowerUBuffer(nCellsX_, 0);
    std::vector<double> receiveLowerVBuffer(nCellsX_, 0);
    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerUBuffer(nCellsX_, 0);
        std::vector<double> sendLowerVBuffer(nCellsX_, 0);
       
        for (int i = 1; i < nCellsX_+1; i++)
        {
            sendLowerUBuffer.at(i-1) = (*discretization_).u(i,1);
            sendLowerVBuffer.at(i-1) = (*discretization_).v(i,1);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendLowerUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendLowerVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveLowerRequests.emplace_back();
        MPI_Irecv(receiveLowerUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLowerRequests.back());

        receiveLowerRequests.emplace_back();
        MPI_Irecv(receiveLowerVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveLowerRequests.back());
    }

    std::vector<double> receiveUpperUBuffer(nCellsX_, 0);
    std::vector<double> receiveUpperVBuffer(nCellsX_, 0);
    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendUpperUBuffer(nCellsX_, 0);
        std::vector<double> sendUpperVBuffer(nCellsX_, 0);
       
        for (int i = 1; i < nCellsX_+1; i++)
        {
            sendUpperUBuffer.at(i-1) = (*discretization_).u(i,nCellsY_);
            sendUpperVBuffer.at(i-1) = (*discretization_).v(i,nCellsY_);
        }
        
        sendRequests.emplace_back();
        MPI_Isend(sendUpperUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        sendRequests.emplace_back();
        MPI_Isend(sendUpperVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 1, MPI_COMM_WORLD, &sendRequests.back());

        receiveUpperRequests.emplace_back();
        MPI_Irecv(receiveUpperUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveUpperRequests.back());

        receiveUpperRequests.emplace_back();
        MPI_Irecv(receiveUpperVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 1, MPI_COMM_WORLD, &receiveUpperRequests.back());
    }


    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Waitall(receiveLeftRequests.size(), receiveLeftRequests.data(), MPI_STATUSES_IGNORE);
        for (int j = 1; j < nCellsY_+1; j++)
        {
            (*discretization_).u(0,j) = receiveLeftUBuffer.at(j-1);
            (*discretization_).v(0,j) = receiveLeftVBuffer.at(j-1);
        }  
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Waitall(receiveRightRequests.size(), receiveRightRequests.data(), MPI_STATUSES_IGNORE);
        for (int j = 1; j < nCellsY_+1; j++)
        {
            (*discretization_).u(nCellsX_+1,j) = receiveRightUBuffer.at(j-1);
            (*discretization_).v(nCellsX_+1,j) = receiveRightVBuffer.at(j-1);
        }  
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Waitall(receiveLowerRequests.size(), receiveLowerRequests.data(), MPI_STATUSES_IGNORE);
        for (int i = 1; i < nCellsX_+1; i++)
        {
            (*discretization_).u(i,0) = receiveLowerUBuffer.at(i-1);
            (*discretization_).v(i,0) = receiveLowerVBuffer.at(i-1);
        }  
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Waitall(receiveUpperRequests.size(), receiveUpperRequests.data(), MPI_STATUSES_IGNORE);
        for (int i = 1; i < nCellsX_+1; i++)
        {
            (*discretization_).u(i,nCellsY_+1) = receiveUpperUBuffer.at(i-1);
            (*discretization_).v(i,nCellsY_+1) = receiveUpperVBuffer.at(i-1);
        }  
    }


    std::vector<MPI_Request> sendDiagonalRequests;
    MPI_Request receiveLeftUpperDiagonalRequest;
    MPI_Request receiveRightLowerDiagonalRequest;

    std::vector<double> receiveLeftUpperDiagonalUBuffer(1, 0);
    std::vector<double> receiveRightLowerDiagonalVBuffer(1, 0);

    if (!partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendLeftUpperDiagonalVBuffer = {(*discretization_).v(1,nCellsY_)};
        
        sendDiagonalRequests.emplace_back();
        MPI_Isend(sendLeftUpperDiagonalVBuffer.data(), 1, MPI_DOUBLE, partitioning_.topNeighbourRankNo()-1, 0, MPI_COMM_WORLD, &sendDiagonalRequests.back());

        MPI_Irecv(receiveLeftUpperDiagonalUBuffer.data(), 1, MPI_DOUBLE, partitioning_.topNeighbourRankNo()-1, 0, MPI_COMM_WORLD, &receiveLeftUpperDiagonalRequest);
    }
    
    if (!partitioning_.ownPartitionContainsRightBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerRightDiagonalUBuffer = {(*discretization_).u(nCellsX_,1)};
        
        sendDiagonalRequests.emplace_back();
        MPI_Isend(sendLowerRightDiagonalUBuffer.data(), 1, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo()+1, 0, MPI_COMM_WORLD, &sendDiagonalRequests.back());

        MPI_Irecv(receiveRightLowerDiagonalVBuffer.data(), 1, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo()+1, 0, MPI_COMM_WORLD, &receiveRightLowerDiagonalRequest);
    }


    MPI_Waitall(sendDiagonalRequests.size(), sendDiagonalRequests.data(), MPI_STATUSES_IGNORE);

    if (!partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&receiveLeftUpperDiagonalRequest, MPI_STATUS_IGNORE);

        (*discretization_).u(0,nCellsY_+1) = receiveLeftUpperDiagonalUBuffer.at(0);
    }

    if (!partitioning_.ownPartitionContainsRightBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&receiveRightLowerDiagonalRequest, MPI_STATUS_IGNORE);

        (*discretization_).v(nCellsX_+1,0) = receiveRightLowerDiagonalVBuffer.at(0);
    }
}

void ComputationParallel::computeTimeStepWidthParallel()
{
    computeTimeStepWidth();

    const double sendDt = dt_;
    double receiveDt;
    
    MPI_Allreduce(&sendDt, &receiveDt, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    dt_ = receiveDt;
}


void ComputationParallel::computePreliminaryVelocities()
{
    int shiftIEndF = 0;
    int shiftJEndG = 0;
    if (!partitioning_.ownPartitionContainsRightBoundary())
        shiftIEndF = 1;
    if (!partitioning_.ownPartitionContainsTopBoundary())
        shiftJEndG = 1;

    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd()+shiftIEndF; i++)
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
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd()+shiftJEndG; j++)
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


void ComputationParallel::receiveAndSendPreliminaryVelocitiesFromAndToOtherProcesses()
{
    std::vector<MPI_Request> sendRequests;
    MPI_Request receiveLeftRequest;
    MPI_Request receiveLowerRequest;

    std::vector<double> receiveLeftFBuffer(nCellsY_, 0);
    std::vector<double> receiveLowerGBuffer(nCellsX_, 0);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Irecv(receiveLeftFBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLeftRequest);
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Irecv(receiveLowerGBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLowerRequest);
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        std::vector<double> sendRightFBuffer(nCellsY_, 0);
       
        for (int j = 1; j < nCellsY_+1; j++)
            sendRightFBuffer.at(j-1) = (*discretization_).f(nCellsX_,j);
        
        sendRequests.emplace_back();
        MPI_Isend(sendRightFBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());
    }

     if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendUpperGBuffer(nCellsX_, 0);
       
        for (int i = 1; i < nCellsX_+1; i++)
            sendUpperGBuffer.at(i-1) = (*discretization_).g(i,nCellsY_);
        
        sendRequests.emplace_back();
        MPI_Isend(sendUpperGBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());
    }


    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Wait(&receiveLeftRequest, MPI_STATUS_IGNORE);

        for (int j = 1; j < nCellsY_+1; j++)
            (*discretization_).f(0,j) = receiveLeftFBuffer.at(j-1);
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&receiveLowerRequest, MPI_STATUS_IGNORE);

        for (int i = 1; i < nCellsX_+1; i++)
            (*discretization_).g(i,0) = receiveLowerGBuffer.at(i-1);
    }
}

void ComputationParallel::computePressure()
{
    (*pressureSolverParallel_).solve();
}

void ComputationParallel::computeVelocities()
{
    int shiftIEndU = 0;
    int shiftJEndV = 0;
    if (!partitioning_.ownPartitionContainsRightBoundary())
        shiftIEndU = 1;
    if (!partitioning_.ownPartitionContainsTopBoundary())
        shiftJEndV = 1;

    for (int i=(*discretization_).uIBegin(); i < (*discretization_).uIEnd()+shiftIEndU; i++)
    {
        for (int j=(*discretization_).uJBegin(); j < (*discretization_).uJEnd(); j++)
        {
            (*discretization_).u(i,j) = (*discretization_).f(i,j) - (dt_/dx_)
                                            * ((*discretization_).p(i+1,j) - (*discretization_).p(i,j));
        }
    }
    for (int i=(*discretization_).vIBegin(); i < (*discretization_).vIEnd(); i++)
    {
        for (int j=(*discretization_).vJBegin(); j < (*discretization_).vJEnd()+shiftJEndV; j++)
        {
            (*discretization_).v(i,j) = (*discretization_).g(i,j) - (dt_/dy_)
                                            * ((*discretization_).p(i,j+1) - (*discretization_).p(i,j));
        }
    }
}
