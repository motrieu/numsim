#include "pressureSolverParallel.h"

#include <mpi.h>


PressureSolverParallel::PressureSolverParallel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, Partitioning partitioning) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning),
    pIBegin_((*discretization_).pIBegin()), pIEnd_((*discretization_).pIEnd()), pJBegin_((*discretization_).pJBegin()), pJEnd_((*discretization_).pJEnd()),
    nCellsX_((*discretization_).nCells()[0]), nCellsY_((*discretization_).nCells()[1])
{
    std::array<int,2> nodeOffset = partitioning.nodeOffset();
    leftAndLowerOffset_ = {(nodeOffset[0]%2) != (nodeOffset[1]%2), (nodeOffset[0]%2) == (nodeOffset[1]%2)};
    rightOffset_ = {(1+leftAndLowerOffset_[0]+nCellsX_)%2, (1+leftAndLowerOffset_[1]+nCellsX_)%2};
    upperOffset_ = {(1+leftAndLowerOffset_[0]+nCellsY_)%2, (1+leftAndLowerOffset_[1]+nCellsY_)%2};

    sendLeftBufferLength_ = {(nCellsY_+1)/2 - leftAndLowerOffset_[0]*(nCellsY_%2), (nCellsY_+1)/2 - leftAndLowerOffset_[1]*(nCellsY_%2)};
    sendLowerBufferLength_ = {(nCellsX_+1)/2 - leftAndLowerOffset_[0]*(nCellsX_%2), (nCellsX_+1)/2 - leftAndLowerOffset_[1]*(nCellsX_%2)};
    sendRightBufferLength_ = {(nCellsY_+1)/2 - rightOffset_[0]*(nCellsY_%2), (nCellsY_+1)/2 - rightOffset_[1]*(nCellsY_%2)};
    sendUpperBufferLength_ = {(nCellsX_+1)/2 - upperOffset_[0]*(nCellsX_%2), (nCellsX_+1)/2 - upperOffset_[1]*(nCellsX_%2)};

    receiveLeftBufferLength_ = {nCellsY_ - sendLeftBufferLength_[0], nCellsY_ - sendLeftBufferLength_[1]};
    receiveLowerBufferLength_ = {nCellsX_ - sendLowerBufferLength_[0], nCellsX_ - sendLowerBufferLength_[1]};
    receiveRightBufferLength_ = {nCellsY_ - sendRightBufferLength_[0], nCellsY_ - sendRightBufferLength_[1]};
    receiveUpperBufferLength_ = {nCellsX_ - sendUpperBufferLength_[0], nCellsX_ - sendUpperBufferLength_[1]};

    numberOfValuesGlobal_ = partitioning_.nCellsGlobal()[0] * partitioning_.nCellsGlobal()[1];
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

void PressureSolverParallel::receiveAndSendPressuresFromAndToOtherProcesses(bool secondHalfStep)
{   
    std::vector<MPI_Request> sendRequests;
    MPI_Request receiveLeftRequest;
    MPI_Request receiveRightRequest;
    MPI_Request receiveLowerRequest;
    MPI_Request receiveUpperRequest;

    std::vector<double> receiveLeftPBuffer(receiveLeftBufferLength_[secondHalfStep], 0);
    std::vector<double> receiveLowerPBuffer(receiveLowerBufferLength_[secondHalfStep], 0);
    std::vector<double> receiveRightPBuffer(receiveRightBufferLength_[secondHalfStep], 0);
    std::vector<double> receiveUpperPBuffer(receiveUpperBufferLength_[secondHalfStep], 0);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        std::vector<double> sendLeftPBuffer(sendLeftBufferLength_[secondHalfStep], 0);
        int k = 0;
        for (int j = pJBegin_+leftAndLowerOffset_[secondHalfStep]; j < pJEnd_; j+=2)
        {
            sendLeftPBuffer[k] = (*discretization_).p(pIBegin_,j);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendLeftPBuffer.data(), sendLeftBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveLeftPBuffer.data(), receiveLeftBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLeftRequest); 
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        std::vector<double> sendLowerPBuffer(sendLowerBufferLength_[secondHalfStep], 0);
        int k = 0;
        for (int i = pIBegin_+leftAndLowerOffset_[secondHalfStep]; i < pIEnd_; i+=2)
        {
            sendLowerPBuffer[k] = (*discretization_).p(i,pJBegin_);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendLowerPBuffer.data(), sendLowerBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveLowerPBuffer.data(), receiveLowerBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveLowerRequest); 
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        std::vector<double> sendRightPBuffer(sendRightBufferLength_[secondHalfStep], 0);
        int k = 0;
        for (int j = pJBegin_+rightOffset_[secondHalfStep]; j < pJEnd_; j+=2)
        {
            sendRightPBuffer[k] = (*discretization_).p(pIEnd_-1,j);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendRightPBuffer.data(), sendRightBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveRightPBuffer.data(), receiveRightBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveRightRequest); 
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        std::vector<double> sendUpperPBuffer(sendUpperBufferLength_[secondHalfStep], 0);
        int k = 0;
        for (int i = pIBegin_+upperOffset_[secondHalfStep]; i < pIEnd_; i+=2)
        {
            sendUpperPBuffer[k] = (*discretization_).p(i,pJEnd_-1);
            k++;
        }
        sendRequests.emplace_back();
        MPI_Isend(sendUpperPBuffer.data(), sendUpperBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &sendRequests.back());

        MPI_Irecv(receiveUpperPBuffer.data(), receiveUpperBufferLength_[secondHalfStep], MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &receiveUpperRequest); 
    }


    MPI_Waitall(sendRequests.size(), sendRequests.data(), MPI_STATUSES_IGNORE);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Wait(&receiveLeftRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int j = pJBegin_+!leftAndLowerOffset_[secondHalfStep]; j < pJEnd_; j+=2)
        {
            (*discretization_).p(pIBegin_-1,j) = receiveLeftPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&receiveLowerRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int i = pIBegin_+!leftAndLowerOffset_[secondHalfStep]; i < pIEnd_; i+=2)
        {
            (*discretization_).p(i,pJBegin_-1) = receiveLowerPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Wait(&receiveRightRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int j = pJBegin_+!rightOffset_[secondHalfStep]; j < pJEnd_; j+=2)
        {
            (*discretization_).p(pIEnd_,j) = receiveRightPBuffer[k];
            k++;
        }
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&receiveUpperRequest, MPI_STATUS_IGNORE);
        int k = 0;
        for (int i = pIBegin_+!upperOffset_[secondHalfStep]; i < pIEnd_; i+=2)
        {
            (*discretization_).p(i,pJEnd_) = receiveUpperPBuffer[k];
            k++;
        }
    }
}

const double PressureSolverParallel::calcResNormSquaredParallel() const
{
    const double sendResNormSquared = calcResNormSquared();
    double receiveResNormSquared;

    MPI_Allreduce(&sendResNormSquared, &receiveResNormSquared, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    return receiveResNormSquared;
}