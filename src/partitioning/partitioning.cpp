#include "partitioning.h"

#include <mpi.h>
#include <iostream>

void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    //calculates partitioning
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks_);

    double globalRatio = static_cast<double>(nCellsGlobal[0]) / static_cast<double>(nCellsGlobal[1]);
    int numPartitionsX;
    int numPartitionsY;

    if (globalRatio >= 1)
    {
        numPartitionsY = std::max(1, roundToInt(std::sqrt(static_cast<double>(numberOfRanks_) / globalRatio)));
        numPartitionsX = numberOfRanks_ / numPartitionsY;
    }
    else
    {
        numPartitionsX = std::max(1, roundToInt(std::sqrt(static_cast<double>(numberOfRanks_) * globalRatio)));
        numPartitionsY = numberOfRanks_ / numPartitionsX;
    }

    assert(numPartitionsY * numPartitionsX == numberOfRanks_);

    int stdNumCellsPerPartitionX = roundToInt(static_cast<double>(nCellsGlobal[0]) / static_cast<double>(numPartitionsX));
    int stdNumCellsPerPartitionY = roundToInt(static_cast<double>(nCellsGlobal[1]) / static_cast<double>(numPartitionsY));

    std::vector<int> numCellsPerPartitionX(numPartitionsX, stdNumCellsPerPartitionX);
    std::vector<int> numCellsPerPartitionY(numPartitionsY, stdNumCellsPerPartitionY);

    int numCellsXOffset = numPartitionsX * stdNumCellsPerPartitionX - nCellsGlobal[0];
    if (numCellsXOffset > 0)
    {
        for (int i=1; i <= numCellsXOffset; i++)
            numCellsPerPartitionX[numCellsPerPartitionX.size()-i] -= 1;
    }
    else if (numCellsXOffset < 0)
    {
        for (int i=1; i <= (-1)*numCellsXOffset; i++)
            numCellsPerPartitionX[numCellsPerPartitionX.size()-i] += 1;
    }
    int numCellsYOffset = numPartitionsY * stdNumCellsPerPartitionY - nCellsGlobal[1];
    if (numCellsYOffset > 0)
    {
        for (int i=1; i <= numCellsYOffset; i++)
            numCellsPerPartitionY[numCellsPerPartitionY.size()-i] -= 1;
    }
    else if (numCellsYOffset < 0)
    {
        for (int i=1; i <= (-1)*numCellsYOffset; i++)
            numCellsPerPartitionY[numCellsPerPartitionY.size()-i] += 1;
    }

    //assigning private variables
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);

    if (ownRank_ < numPartitionsX)
        rankLower_ = -1;
    else
        rankLower_ = ownRank_ - numPartitionsX;
    if (ownRank_ >= ((numPartitionsY-1) * numPartitionsX))
        rankUpper_ = -1;
    else
        rankUpper_ = ownRank_ + numPartitionsX;
    if ((ownRank_ % numPartitionsX) == 0)
        rankLeft_ = -1;
    else
        rankLeft_ = ownRank_ - 1;
    if ((ownRank_ % numPartitionsX) == numPartitionsX - 1)
        rankRight_ = -1;
    else
        rankRight_ = ownRank_ + 1;

    ownPartitionContainsBottomBoundary_ = rankLower_ == -1;
    ownPartitionContainsTopBoundary_ = rankUpper_ == -1;
    ownPartitionContainsLeftBoundary_ = rankLeft_ == -1;
    ownPartitionContainsRightBoundary_ = rankRight_ == -1;

    int partitionIndexX = ownRank_%numPartitionsX;
    int partitionIndexY = ownRank_/numPartitionsX;
    nCellsLocal_ = {numCellsPerPartitionX[partitionIndexX], numCellsPerPartitionY[partitionIndexY]};
    nCellsGlobal_ = nCellsGlobal;

    int nodeOffsetX = 0;
    int nodeOffsetY = 0;
    for (int i = 0; i < partitionIndexX; i++)
        nodeOffsetX += numCellsPerPartitionX[i];
    for (int j = 0; j < partitionIndexY; j++)
        nodeOffsetY += numCellsPerPartitionY[j];
    nodeOffset_ = {nodeOffsetX, nodeOffsetY};
}

std::array<int, 2> Partitioning::nCellsLocal() const
{
    return nCellsLocal_;
}

std::array<int, 2> Partitioning::nCellsGlobal() const
{
    return nCellsGlobal_;
}

int Partitioning::ownRankNo() const
{
    return ownRank_;
}

int Partitioning::nRanks() const
{
    return numberOfRanks_;
}

bool Partitioning::ownPartitionContainsBottomBoundary() const
{
    return ownPartitionContainsBottomBoundary_;
}

bool Partitioning::ownPartitionContainsTopBoundary() const
{
    return ownPartitionContainsTopBoundary_;
}

bool Partitioning::ownPartitionContainsLeftBoundary() const
{
    return ownPartitionContainsLeftBoundary_;
}

bool Partitioning::ownPartitionContainsRightBoundary() const
{
    return ownPartitionContainsRightBoundary_;
}

int Partitioning::leftNeighbourRankNo() const
{
    return rankLeft_;
}

int Partitioning::rightNeighbourRankNo() const
{
    return rankRight_;
}

int Partitioning::topNeighbourRankNo() const
{
    return rankUpper_;
}

int Partitioning::bottomNeighbourRankNo() const
{
    return rankLower_;
}

std::array<int, 2> Partitioning::nodeOffset() const
{
    return nodeOffset_;
}

int Partitioning::roundToInt(double x)
{
    int xInt = x;
    if ((x - std::floor(x)) >= 0.5)
        xInt += 1;
    return xInt;
}