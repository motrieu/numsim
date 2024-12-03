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
    MPI_Request leftRequest;
    MPI_Request rightRequest;
    MPI_Request lowerRequest;
    MPI_Request upperRequest;

    std::vector<double> leftPBuffer(nCellsY_); //receiveLeftBufferLength_[secondHalfStep], 0);
    std::vector<double> lowerPBuffer(nCellsX_); //receiveLowerBufferLength_[secondHalfStep], 0);
    std::vector<double> rightPBuffer(nCellsY_); //receiveRightBufferLength_[secondHalfStep], 0);
    std::vector<double> upperPBuffer(nCellsX_); //receiveUpperBufferLength_[secondHalfStep], 0);

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        //int k = 0;
        for (int j = pJBegin_/*+leftAndLowerOffset_[secondHalfStep]*/; j < pJEnd_; j++)//+=2)
        {
            leftPBuffer.at(j-pJBegin_) = (*discretization_).p(pIBegin_,j);
            //k++;
        }

        MPI_Isend(leftPBuffer.data(), nCellsY_/*sendLeftBufferLength_[secondHalfStep]*/, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftRequest);

        MPI_Irecv(leftPBuffer.data(), nCellsY_/*receiveLeftBufferLength_[secondHalfStep]*/, MPI_DOUBLE, partitioning_.leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &leftRequest); 
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        for (int i = pIBegin_; i < pIEnd_; i++)
        {
            lowerPBuffer.at(i-pIBegin_) = (*discretization_).p(i,pJBegin_);
        }
        
        MPI_Isend(lowerPBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerRequest);

        MPI_Irecv(lowerPBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &lowerRequest); 
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        for (int j = pJBegin_; j < pJEnd_; j++)
        {
            rightPBuffer.at(j-pJBegin_) = (*discretization_).p(pIEnd_-1,j);
        }
        
        MPI_Isend(rightPBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightRequest);

        MPI_Irecv(rightPBuffer.data(), nCellsY_, MPI_DOUBLE, partitioning_.rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &rightRequest); 
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        for (int i = pIBegin_; i < pIEnd_; i++)
        {
            upperPBuffer.at(i-pIBegin_) = (*discretization_).p(i,pJEnd_-1);
        }
        
        MPI_Isend(upperPBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperRequest);

        MPI_Irecv(upperPBuffer.data(), nCellsX_, MPI_DOUBLE, partitioning_.topNeighbourRankNo(), 0, MPI_COMM_WORLD, &upperRequest); 
    }

    if (!partitioning_.ownPartitionContainsLeftBoundary())
    {
        MPI_Wait(&leftRequest, MPI_STATUS_IGNORE);
        for (int j = pJBegin_; j < pJEnd_; j++)
        {
            (*discretization_).p(pIBegin_-1,j) = leftPBuffer.at(j-pJBegin_);
        }
    }

    if (!partitioning_.ownPartitionContainsBottomBoundary())
    {
        MPI_Wait(&lowerRequest, MPI_STATUS_IGNORE);
        for (int i = pIBegin_; i < pIEnd_; i++)
        {
            (*discretization_).p(i,pJBegin_-1) = lowerPBuffer.at(i-pIBegin_);
        }
    }

    if (!partitioning_.ownPartitionContainsRightBoundary())
    {
        MPI_Wait(&rightRequest, MPI_STATUS_IGNORE);
        for (int j = pJBegin_; j < pJEnd_; j++)
        {
            (*discretization_).p(pIEnd_,j) = rightPBuffer.at(j-pJBegin_);
        }
    }

    if (!partitioning_.ownPartitionContainsTopBoundary())
    {
        MPI_Wait(&upperRequest, MPI_STATUS_IGNORE);
        for (int i = pIBegin_; i < pIEnd_; i++)
        {
            (*discretization_).p(i,pJEnd_) = upperPBuffer.at(i-pIBegin_);
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