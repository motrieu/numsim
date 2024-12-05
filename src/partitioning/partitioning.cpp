#include "partitioning.h"

#include <mpi.h>
#include <iostream>

void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    MPI_Comm_size(MPI_COMM_WORLD, &numberOfRanks_);

    // calculate ratio describing the shape of domain,
    // is 1 for same discretization in x and y (quadratic) and !=1 for rectangular discretizations
    // the aim is to fit the number of partitions in x and y direction to the number of cells in x and y direction and therefore to the global ratio
    double globalRatio = static_cast<double>(nCellsGlobal[0]) / static_cast<double>(nCellsGlobal[1]);
    int numPartitionsX;
    int numPartitionsY;

    // If the discretization is more expanded in horizontal direction (horizontal rectangle), the priority is to
    // have more partitions in x direction. Therefore, a first guess for the number of partitions in y direction is
    // calculated which might then be reduced till a valid decomposition has been found.
    // The first guess is based on two equations:
    //      1. as demanded above, the ratio of partitions in x and y should be similar to the ratio of cells in x and y:
    //              nCellsX/nCellsY = nPartX/nPartY
    //      2. number of partitions should be equivalent to the number of MPI ranks:
    //              nPartX*nPartY = nRanks
    //     ==> nPartY = sqrt(nRanks*nCellsY/nCellsX)
    // Note that the second equation is not guaranteed to be fulfilled due to rounding of the first guess.
    // To ensure the second equation the non-priority number of partitions in y is reduced until it is fulfilled.
    if (globalRatio >= 1)
    {
        numPartitionsY = std::max(1, roundToInt(std::sqrt(static_cast<double>(numberOfRanks_) / globalRatio)));
        numPartitionsX = numberOfRanks_ / numPartitionsY;
        while (numPartitionsX*numPartitionsY != numberOfRanks_)
        {
            numPartitionsY -= 1;
            numPartitionsX = numberOfRanks_ / numPartitionsY;
        }
    }
    // If the discretization is more expanded in vertical direction (vertical rectangle), the priority is to
    // have more partitions in y direction. Therefore, a first guess for the number of partitions in x direction is
    // calculated which might then be reduced till a valid decomposition has been found.
    // The same equations as above yield:
    //      nPartX = sqrt(nRanks*nCellsX/nCellsY)
    else
    {
        numPartitionsX = std::max(1, roundToInt(std::sqrt(static_cast<double>(numberOfRanks_) * globalRatio)));
        numPartitionsY = numberOfRanks_ / numPartitionsX;
        while (numPartitionsX*numPartitionsY != numberOfRanks_)
        {
            numPartitionsX -= 1;
            numPartitionsY = numberOfRanks_ / numPartitionsX;
        }
    }

    assert(numPartitionsY * numPartitionsX == numberOfRanks_);

    // calculates standard number of cells per partition in x and y
    int stdNumCellsPerPartitionX = roundToInt(static_cast<double>(nCellsGlobal[0]) / static_cast<double>(numPartitionsX));
    int stdNumCellsPerPartitionY = roundToInt(static_cast<double>(nCellsGlobal[1]) / static_cast<double>(numPartitionsY));

    // it might be that a multiple of the standard value does not match the global number of cells
    // if using the standard value would lead to too many global cells, the number of cells is reduced by one each in
    // the appropriate number of last partitions
    // if using the standard value would lead to too few global cells, the number of cells is increased by one each in
    // the appropriate number of last partitions
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

    // assigning own rank and neighbouring rank numbers

    MPI_Comm_rank(MPI_COMM_WORLD, &ownRank_);

    // if the ownRank is in the lowest row, it has no lower neighbour, i.e. the number of the lower rank is set to -1
    if (ownRank_ < numPartitionsX)
        rankLower_ = -1;
    // otherwise the number of the lower rank is calculated by reducing the own rank number by one whole row
    else
        rankLower_ = ownRank_ - numPartitionsX;

    // analogously to above for the upper, left and right neigbours 
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

    // position of partition in global context in rank indexing
    int partitionIndexX = ownRank_%numPartitionsX;
    int partitionIndexY = ownRank_/numPartitionsX;

    nCellsLocal_ = {numCellsPerPartitionX[partitionIndexX], numCellsPerPartitionY[partitionIndexY]};
    nCellsGlobal_ = nCellsGlobal;

    // position of partition in global context in cell indexing
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