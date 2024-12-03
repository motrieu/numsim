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
            //receiveAndSendDiagonalPressureFromAndToOtherProcess();

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
    //receiveAndSendDiagonalPressureFromAndToOtherProcess();
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
    MPI_Request diagonalRequest;

    if ((!partitioning_.ownPartitionContainsLeftBoundary()) && (!partitioning_.ownPartitionContainsBottomBoundary()))
    {
        std::vector<double> sendLeftLowerDiagonalPBuffer = {(*discretization_).p(1,1)};

        MPI_Isend(sendLeftLowerDiagonalPBuffer.data(), 1, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo()-1, 0, MPI_COMM_WORLD, &diagonalRequest);
    }

    std::vector<double> receiveRightUpperDiagonalPBuffer(1);
    if ((!partitioning_.ownPartitionContainsRightBoundary()) && (!partitioning_.ownPartitionContainsTopBoundary()))
    {
        MPI_Irecv(receiveRightUpperDiagonalPBuffer.data(), 1, MPI_DOUBLE, partitioning_.topNeighbourRankNo()+1, 0, MPI_COMM_WORLD, &diagonalRequest);
    }


    if ((!partitioning_.ownPartitionContainsRightBoundary()) && (!partitioning_.ownPartitionContainsTopBoundary()))
    {
        MPI_Wait(&diagonalRequest, MPI_STATUS_IGNORE);

        (*discretization_).p(nCellsX_+1,nCellsY_+1) = receiveRightUpperDiagonalPBuffer[0];
    }
}

void ComputationParallel::receiveAndSendVelocitiesFromAndToOtherProcesses()
{
    MPI_Request leftURequest;
    MPI_Request leftVRequest;
    MPI_Request rightURequest;
    MPI_Request rightVRequest;
    MPI_Request lowerURequest;
    MPI_Request lowerVRequest;
    MPI_Request upperURequest;
    MPI_Request upperVRequest;

    std::vector<double> leftUBuffer(nCellsY_);
    std::vector<double> leftVBuffer(nCellsY_);
    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        for (int j = 1; j < nCellsY_+1; j++)
        {
            leftUBuffer.at(j-1) = (*discretization_).u(1,j);
            leftVBuffer.at(j-1) = (*discretization_).v(1,j);
        }
        
        MPI_Isend(leftUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftURequest);
        MPI_Isend(leftVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftVRequest);

        MPI_Irecv(leftUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftURequest);
        MPI_Irecv(leftVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftVRequest);
    }

    std::vector<double> rightUBuffer(nCellsY_);
    std::vector<double> rightVBuffer(nCellsY_);
    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        for (int j = 1; j < nCellsY_+1; j++)
        {
            rightUBuffer.at(j-1) = (*discretization_).u(nCellsX_,j);
            rightVBuffer.at(j-1) = (*discretization_).v(nCellsX_,j);
        }
        
        MPI_Isend(rightUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightURequest);
        MPI_Isend(rightVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightVRequest);

        MPI_Irecv(rightUBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightURequest);
        MPI_Irecv(rightVBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightVRequest);
    }

    std::vector<double> lowerUBuffer(nCellsX_);
    std::vector<double> lowerVBuffer(nCellsX_);
    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = 1; i < nCellsX_+1; i++)
        {
            lowerUBuffer.at(i-1) = (*discretization_).u(i,1);
            lowerVBuffer.at(i-1) = (*discretization_).v(i,1);
        }
        
        MPI_Isend(lowerUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerURequest);
        MPI_Isend(lowerVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerVRequest);

        MPI_Irecv(lowerUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerURequest);
        MPI_Irecv(lowerVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerVRequest);
    }

    std::vector<double> upperUBuffer(nCellsX_);
    std::vector<double> upperVBuffer(nCellsX_);
    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = 1; i < nCellsX_+1; i++)
        {
            upperUBuffer.at(i-1) = (*discretization_).u(i,nCellsY_);
            upperVBuffer.at(i-1) = (*discretization_).v(i,nCellsY_);
        }
        
        MPI_Isend(upperUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperURequest);
        MPI_Isend(upperVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperVRequest);

        MPI_Irecv(upperUBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperURequest);
        MPI_Irecv(upperVBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperVRequest);
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Wait(&leftURequest, MPI_STATUSES_IGNORE);
        MPI_Wait(&leftVRequest, MPI_STATUSES_IGNORE);
        for (int j = 1; j < nCellsY_+1; j++)
        {
            (*discretization_).u(0,j) = leftUBuffer.at(j-1);
            (*discretization_).v(0,j) = leftVBuffer.at(j-1);
        }  
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Wait(&rightURequest, MPI_STATUSES_IGNORE);
        MPI_Wait(&rightVRequest, MPI_STATUSES_IGNORE);
        for (int j = 1; j < nCellsY_+1; j++)
        {
            (*discretization_).u(nCellsX_+1,j) = rightUBuffer.at(j-1);
            (*discretization_).v(nCellsX_+1,j) = rightVBuffer.at(j-1);
        }  
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&lowerURequest, MPI_STATUSES_IGNORE);
        MPI_Wait(&lowerVRequest, MPI_STATUSES_IGNORE);
        for (int i = 1; i < nCellsX_+1; i++)
        {
            (*discretization_).u(i,0) = lowerUBuffer.at(i-1);
            (*discretization_).v(i,0) = lowerVBuffer.at(i-1);
        }  
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&upperURequest, MPI_STATUSES_IGNORE);
        MPI_Wait(&upperVRequest, MPI_STATUSES_IGNORE);
        for (int i = 1; i < nCellsX_+1; i++)
        {
            (*discretization_).u(i,nCellsY_+1) = upperUBuffer.at(i-1);
            (*discretization_).v(i,nCellsY_+1) = upperVBuffer.at(i-1);
        }  
    }


    /*MPI_Request leftUpperDiagonalRequest;
    MPI_Request rightLowerDiagonalRequest;

    std::vector<double> receiveLeftUpperDiagonalUBuffer(1);
    std::vector<double> receiveRightLowerDiagonalVBuffer(1);

    if (!partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendLeftUpperDiagonalVBuffer = {(*discretization_).v(1,nCellsY_)};
        
        MPI_Isend(sendLeftUpperDiagonalVBuffer.data(), 1, MPI_DOUBLE, partitioning_.topNeighbourRankNo()-1, 0, MPI_COMM_WORLD, &leftUpperDiagonalRequest);

        MPI_Irecv(receiveLeftUpperDiagonalUBuffer.data(), 1, MPI_DOUBLE, partitioning_.topNeighbourRankNo()-1, 0, MPI_COMM_WORLD, &leftUpperDiagonalRequest);
    }
    
    if (!partitioning_.ownPartitionContainsRightBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerRightDiagonalUBuffer = {(*discretization_).u(nCellsX_,1)};
        
        MPI_Isend(sendLowerRightDiagonalUBuffer.data(), 1, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo()+1, 0, MPI_COMM_WORLD, &rightLowerDiagonalRequest);

        MPI_Irecv(receiveRightLowerDiagonalVBuffer.data(), 1, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo()+1, 0, MPI_COMM_WORLD, &rightLowerDiagonalRequest);
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary() && !partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&leftUpperDiagonalRequest, MPI_STATUS_IGNORE);

        (*discretization_).u(0,nCellsY_+1) = receiveLeftUpperDiagonalUBuffer[0];
    }

    if (!partitioning_.ownPartitionContainsRightBoundary() && !partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&rightLowerDiagonalRequest, MPI_STATUS_IGNORE);

        (*discretization_).v(nCellsX_+1,0) = receiveRightLowerDiagonalVBuffer[0];
    }*/
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
    MPI_Request leftRequest;
    MPI_Request rightRequest;
    MPI_Request lowerRequest;
    MPI_Request upperRequest;

    std::vector<double> leftFBuffer(nCellsY_);
    std::vector<double> lowerGBuffer(nCellsX_);
    std::vector<double> rightFBuffer(nCellsY_);
    std::vector<double> upperGBuffer(nCellsX_);

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        for (int j = 1; j < nCellsY_+1; j++)
            rightFBuffer.at(j-1) = (*discretization_).f(nCellsX_,j);
        
        MPI_Isend(rightFBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightRequest);
    
        MPI_Irecv(rightFBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightRequest);
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = 1; i < nCellsX_+1; i++)
            upperGBuffer.at(i-1) = (*discretization_).g(i,nCellsY_);
        
        MPI_Isend(upperGBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperRequest);

        MPI_Irecv(upperGBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperRequest);
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        for (int j = 1; j < nCellsY_+1; j++)
            leftFBuffer.at(j-1) = (*discretization_).f(0,j);
        
        MPI_Isend(leftFBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftRequest);

        MPI_Irecv(leftFBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftRequest);
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = 1; i < nCellsX_+1; i++)
            lowerGBuffer.at(i-1) = (*discretization_).g(i,0);
        
        MPI_Isend(lowerGBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerRequest);
        
        MPI_Irecv(lowerGBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerRequest);
    }


    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Wait(&rightRequest, MPI_STATUS_IGNORE);

        for (int j = 1; j < nCellsY_+1; j++)
            (*discretization_).f(nCellsX_+1,j) = rightFBuffer.at(j-1);
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&upperRequest, MPI_STATUS_IGNORE);

        for (int i = 1; i < nCellsX_+1; i++)
            (*discretization_).g(i,nCellsY_+1) = upperGBuffer.at(i-1);
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Wait(&leftRequest, MPI_STATUS_IGNORE);

        for (int j = 1; j < nCellsY_+1; j++)
            (*discretization_).f(0,j) = leftFBuffer.at(j-1);
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&lowerRequest, MPI_STATUS_IGNORE);

        for (int i = 1; i < nCellsX_+1; i++)
            (*discretization_).g(i,0) = lowerGBuffer.at(i-1);
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
