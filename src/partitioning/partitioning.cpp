#include "partitioning.h"

void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    int nPartitions;
    MPI_Comm_size(MPI_COMM_WORLD, &nPartitions);

    double globalRatio = static_cast<double>(nCellsGlobal[0]) / static_cast<double>(nCellsGlobal[1]);
    int numberPartitionsX;
    int numberPartitionsY;

    if (globalRatio >= 1)
    {
        numberPartitionsY = std::round(std::sqrt(static_cast<double>(nPartitions) / globalRatio));
        numberPartitionsX = nPartitions / numberPartitionsY;
    }
    else
    {
        numberPartitionsX = std::round(std::sqrt(static_cast<double>(nPartitions) * globalRatio));
        numberPartitionsY = nPartitions / numberPartitionsX;
    }

    assert(numberPartitionsY * numberPartitionsX == nPartitions);

    int stdNumCellsPerPartitionX = std::round(static_cast<double>(nCellsGlobal[0]) / static_cast<double>(numberPartitionsX));
    int stdNumCellsPerPartitionY = std::round(static_cast<double>(nCellsGlobal[1]) / static_cast<double>(numberPartitionsY));

    std::array<int> numCellsPerPartitionX(numberPartitionsX, stdNumCellsPerPartitionX);
    std::array<int> numCellsPerPartitionY(numberPartitionsY, stdNumCellsPerPartitionY);

    if (numberPartitionsX * stdNumCellsPerPartitionX > nCellsGlobal[0])
    {
    }
    else if (numberPartitionsX * stdNumCellsPerPartitionX < nCellsGlobal[0])
}

std::array<int, 2> nCellsLocal() const
{
}

std::array<int, 2> nCellsGlobal() const
{
}

int ownRankNo() const
{
}

int nRanks() const
{
}

bool ownPartitionContainsBottomBoundary() const
{
}

bool ownPartitionContainsTopBoundary() const
{
}

bool ownPartitionContainsLeftBoundary() const
{
}

bool ownPartitionContainsRightBoundary() const
{
}

int leftNeighbourRankNo() const
{
}

int rightNeighbourRankNo() const
{
}

int topNeighbourRankNo() const
{
}

int bottomNeighbourRankNo() const
{
}

std::array<int, 2> nodeOffset() const
{
}