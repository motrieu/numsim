#include "pressureSolverParallel.h"

#include <mpi.h>


PressureSolverParallel::PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning),
    pIBegin_((*discretization_).pIBegin()), pIEnd_((*discretization_).pIEnd()), pJBegin_((*discretization_).pJBegin()), pJEnd_((*discretization_).pJEnd()),
    nCellsX_((*discretization_).nCells()[0]), nCellsY_((*discretization_).nCells()[1])
{
    std::array<int,2> nodeOffset = partitioning.nodeOffset();
    leftAndLowerOffset_ = (nodeOffset[0]%2) != (nodeOffset[1]%2);
    rightOffset_ = (1+leftAndLowerOffset_+nCellsX_)%2;
    upperOffset_ = (1+leftAndLowerOffset_+nCellsY_)%2;
}

void PressureSolverParallel::setBoundaryValuesOnDirichletParallel()
{
    if (partitioning_.ownPartitionContainsLeftBoundary())
    {
        for (int j = pJBegin_; j < pJEnd_; j++)
            (*discretization_).p(pIBegin_-1,j) = (*discretization_).p(pIBegin_,j);
    }

    if (partitioning_.ownPartitionContainsRightBoundary())
    {
        for (int j = pJBegin_; j < pJEnd_; j++)
            (*discretization_).p(pIEnd_,j) = (*discretization_).p(pIEnd_-1,j);
    }

    if (partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = pIBegin_; i < pIEnd_; i++)
            (*discretization_).p(i,pJBegin_-1) = (*discretization_).p(i,pJBegin_);
    }

    if (partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = pIBegin_; i < pIEnd_; i++)
            (*discretization_).p(i,pJEnd_) = (*discretization_).p(i,pJEnd_-1);
    }
}

void PressureSolverParallel::receiveAndSendPressuresFromAndToOtherProcesses(bool leftAndLowerOffset, bool rightOffset, bool upperOffset)
{   
    std::vector<MPI_Request> sendRequests;
    std::vector<MPI_Request> receiveRequests;

    int sendLeftBufferLength = (nCellsY_+1)/2 - leftAndLowerOffset*(nCellsY_%2);
    int sendLowerBufferLength = (nCellsX_+1)/2 - leftAndLowerOffset*(nCellsX_%2);
    int sendRightBufferLength = (nCellsY_+1)/2 - rightOffset*(nCellsY_%2);
    int sendUpperBufferLength = (nCellsX_+1)/2 - upperOffset*(nCellsX_%2);

    int receiveLeftBufferLength = nCellsY_ - sendLeftBufferLength;
    int receiveLowerBufferLength = nCellsX_ - sendLowerBufferLength;
    int receiveRightBufferLength = nCellsY_ - sendRightBufferLength;
    int receiveUpperBufferLength = nCellsX_ - sendUpperBufferLength;

    std::vector<double> receiveLeftPBuffer(receiveLeftBufferLength);
    std::vector<double> receiveLowerPBuffer(receiveLowerBufferLength);
    std::vector<double> receiveRightPBuffer(receiveRightBufferLength);
    std::vector<double> receiveUpperPBuffer(receiveUpperBufferLength);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        std::vector<double> sendLeftPBuffer(sendLeftBufferLength);
        int k = 0;
        for (int j = pJBegin_+leftAndLowerOffset; j < pJEnd_; j+=2)
        {
            sendLeftPBuffer[k] = (*discretization_).p(pIBegin_,j);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendLeftPBuffer.data(), sendLeftBufferLength, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        receiveRequests.emplace_back();
        MPI_Irecv(receiveLeftPBuffer.data(), receiveLeftBufferLength, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRequests.back()); 
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerPBuffer(sendLowerBufferLength);
        int k = 0;
        for (int i = pIBegin_+leftAndLowerOffset; i < pIEnd_; i+=2)
        {
            sendLowerPBuffer[k] = (*discretization_).p(i,pJBegin_);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendLowerPBuffer.data(), sendLowerBufferLength, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        receiveRequests.emplace_back();
        MPI_Irecv(receiveLowerPBuffer.data(), receiveLowerBufferLength, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRequests.back()); 
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        std::vector<double> sendRightPBuffer(sendRightBufferLength);
        int k = 0;
        for (int j = pJBegin_+rightOffset; j < pJEnd_; j+=2)
        {
            sendRightPBuffer[k] = (*discretization_).p(pIEnd_-1,j);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendRightPBuffer.data(), sendRightBufferLength, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        receiveRequests.emplace_back();
        MPI_Irecv(receiveRightPBuffer.data(), receiveRightBufferLength, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRequests.back()); 
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendUpperPBuffer(sendUpperBufferLength);
        int k = 0;
        for (int i = pIBegin_+upperOffset; i < pIEnd_; i+=2)
        {
            sendUpperPBuffer[k] = (*discretization_).p(i,pJEnd_-1);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendUpperPBuffer.data(), sendUpperBufferLength, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        receiveRequests.emplace_back();
        MPI_Irecv(receiveUpperPBuffer.data(), receiveUpperBufferLength, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRequests.back()); 
    }


    MPI_Waitall(receiveRequests.size(), receiveRequests.data(), MPI_STATUSES_IGNORE);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        int k = 0;
        for (int j = pJBegin_+!leftAndLowerOffset; j < pJEnd_; j+=2)
        {
            (*discretization_).p(pIBegin_-1,j) = receiveLeftPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        int k = 0;
        for (int i = pIBegin_+!leftAndLowerOffset; i < pIEnd_; i+=2)
        {
            (*discretization_).p(i,pJBegin_-1) = receiveLowerPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        int k = 0;
        for (int j = pJBegin_+!rightOffset; j < pJEnd_; j+=2)
        {
            (*discretization_).p(pIEnd_,j) = receiveRightPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        int k = 0;
        for (int i = pIBegin_+!upperOffset; i < pIEnd_; i+=2)
        {
            (*discretization_).p(i,pJEnd_) = receiveUpperPBuffer[k];
            k++;
        }
    }

    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);
}

const double PressureSolverParallel::calcResNormSquaredParallel() const
{
    const double sendResNormSquared = calcResNormSquared();
    double receiveResNormSquared;

    MPI_Allreduce(&sendResNormSquared, &receiveResNormSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return receiveResNormSquared;
}