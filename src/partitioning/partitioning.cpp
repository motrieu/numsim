#include "partitioning.h"

#include <iostream>

void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    int nPartitions = 48;
    //MPI_Comm_size(MPI_COMM_WORLD, &nPartitions);

    double globalRatio = static_cast<double>(nCellsGlobal[0]) / static_cast<double>(nCellsGlobal[1]);
    int numPartitionsX;
    int numPartitionsY;

    if (globalRatio >= 1)
    {
        numPartitionsY = std::max(1, roundToInt(std::sqrt(static_cast<double>(nPartitions) / globalRatio)));
        numPartitionsX = nPartitions / numPartitionsY;
    }
    else
    {
        numPartitionsX = std::max(1, roundToInt(std::sqrt(static_cast<double>(nPartitions) * globalRatio)));
        numPartitionsY = nPartitions / numPartitionsX;
    }

    assert(numPartitionsY * numPartitionsX == nPartitions);

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

    for (int i=0; i < numCellsPerPartitionX.size(); i++)
        std::cout << numCellsPerPartitionX[i] << ", ";
    std::cout << std::endl;
    for (int i=0; i < numCellsPerPartitionY.size(); i++)
        std::cout << numCellsPerPartitionY[i] << ", ";
    std::cout << std::endl;
}

/*std::array<int, 2> nCellsLocal() const
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
}*/

int Partitioning::roundToInt(double x)
{
    int xInt = x;
    if ((x - std::floor(x)) >= 0.5)
        xInt += 1;
    return xInt;
}