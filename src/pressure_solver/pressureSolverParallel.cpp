#include "pressureSolverParallel.h"


PressureSolverParallel::PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning)
{
    std::array<int,2> nodeOffset = partitioning.nodeOffset();
    startOffset_ = (nodeOffset[0]%2) != (nodeOffset[1]%2);
}

void PressureSolverParallel::setBoundaryValuesOnDirichletParallel()
{
    int pIBegin = (*discretization_).pIBegin();
    int pIEnd = (*discretization_).pIEnd();
    int pJBegin = (*discretization_).pJBegin();
    int pJEnd = (*discretization_).pJEnd();

    if (partitioning_.ownPartitionContainsLeftBoundary())
    {
        for (int j = pJBegin; j < pJEnd; j++)
            (*discretization_).p(pIBegin-1,j) = (*discretization_).p(pIBegin,j);
    }

    if (partitioning_.ownPartitionContainsRightBoundary())
    {
        for (int j = pJBegin; j < pJEnd; j++)
            (*discretization_).p(pIEnd,j) = (*discretization_).p(pIEnd-1,j);
    }

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = pIBegin; i < pIEnd; i++)
            (*discretization_).p(i,pJBegin-1) = (*discretization_).p(i,pJBegin);
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = pIBegin; i < pIEnd; i++)
            (*discretization_).p(i,pJEnd) = (*discretization_).p(i,pJEnd-1);
    }
}

void PressureSolverParallel::receiveAndSendPressuresFromAndToOtherProcesses()
{   
    std::vector<MPI_Request> sendRequests;
    MPI_Request receiveLeftRequest;
    MPI_Request receiveRightRequest;
    MPI_Request receiveLowerRequest;
    MPI_Request receiveUpperRequest;

    int nCellsX = (*discretization_).nCells()[0];
    int nCellsY = (*discretization_).nCells()[1];

    int pIBegin = (*discretization_).pIBegin();
    int pIEnd = (*discretization_).pIEnd();
    int pJBegin = (*discretization_).pJBegin();
    int pJEnd = (*discretization_).pJEnd();

    int rightOffset = (1+startOffset_+nCellsX)%2;
    int upperOffset = (1+startOffset_+nCellsY)%2;

    int sendLeftBufferLength = (nCellsY+1)/2 - startOffset_*(nCellsY%2);
    int sendLowerBufferLength = (nCellsX+1)/2 - startOffset_*(nCellsX%2);
    int sendRightBufferLength = (nCellsY+1)/2 - rightOffset*(nCellsY%2);
    int sendUpperBufferLength = (nCellsX+1)/2 - upperOffset*(nCellsX%2);

    int receiveLeftBufferLength = nCellsY - sendLeftBufferLength;
    int receiveLowerBufferLength = nCellsX - sendLowerBufferLength;
    int receiveRightBufferLength = nCellsY - sendRightBufferLength;
    int receiveUpperBufferLength = nCellsX - sendUpperBufferLength;

    std::vector<double> receiveLeftPBuffer(receiveLeftBufferLength);
    std::vector<double> receiveLowerPBuffer(receiveLowerBufferLength);
    std::vector<double> receiveRightPBuffer(receiveRightBufferLength);
    std::vector<double> receiveUpperPBuffer(receiveUpperBufferLength);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        std::vector<double> sendLeftPBuffer(sendLeftBufferLength);
        int k = 0;
        for (int j = pJBegin+startOffset_; j < pJEnd; j+=2)
        {
            sendLeftPBuffer[k] = (*discretization_).p(pIBegin,j);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendLeftPBuffer.data(), sendLeftBufferLength, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveLeftPBuffer.data(), receiveLeftBufferLength, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLeftRequest); 
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerPBuffer(sendLowerBufferLength);
        int k = 0;
        for (int i = pIBegin+startOffset_; i < pIEnd; i+=2)
        {
            sendLowerPBuffer[k] = (*discretization_).p(i,pJBegin);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendLowerPBuffer.data(), sendLowerBufferLength, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveLowerPBuffer.data(), receiveLowerBufferLength, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLowerRequest); 
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        std::vector<double> sendRightPBuffer(sendRightBufferLength);
        int k = 0;
        for (int j = pJBegin+rightOffset; j < pJEnd; j+=2)
        {
            sendRightPBuffer[k] = (*discretization_).p(pIEnd-1,j);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendRightPBuffer.data(), sendRightBufferLength, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveRightPBuffer.data(), receiveRightBufferLength, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRightRequest); 
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendUpperPBuffer(sendUpperBufferLength);
        int k = 0;
        for (int i = pIBegin+upperOffset; i < pIEnd; i+=2)
        {
            sendUpperPBuffer[k] = (*discretization_).p(i,pJEnd-1);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendUpperPBuffer.data(), sendUpperBufferLength, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveUpperPBuffer.data(), receiveUpperBufferLength, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveUpperRequest); 
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Wait(&receiveLeftRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int j = pJBegin+!startOffset_; j < pJEnd; j+=2)
        {
            (*discretization_).p(pIBegin-1,j) = receiveLeftPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&receiveLowerRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int i = pIBegin+!startOffset_; i < pIEnd; i+=2)
        {
            (*discretization_).p(i,pJBegin-1) = receiveLowerPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Wait(&receiveRightRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int j = pJBegin+!rightOffset; j < pJEnd; j+=2)
        {
            (*discretization_).p(pJEnd,j) = receiveRightPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&receiveUpperRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int i = pIBegin+!upperOffset; i < pIEnd; i+=2)
        {
            (*discretization_).p(i,pJEnd) = receiveUpperPBuffer[k];
            k++;
        }
    }

    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
}
