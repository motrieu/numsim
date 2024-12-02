#include "computation/computationParallel.h"

#include <array>
#include <vector>
#include <mpi.h>
#include <unistd.h>
#include <iostream>

void debugComm()
{
  int numRanks;
  MPI_Comm_size(MPI_COMM_WORLD, &numRanks);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  int leftMessage = 0;
  int rightMessage = 0;

  int receiveLeftBuffer = 0;
  int receiveRightBuffer = 0;

  int sendLeftBuffer = 1;
  int sendRightBuffer = 1;

  std::vector<MPI_Request> sendRequest;
  std::vector<MPI_Request> receiveRequest;

  if (rank > 0)
  {
    sendRequest.emplace_back();
    receiveRequest.emplace_back();

    MPI_Isend(&sendLeftBuffer, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &sendRequest.back());
    MPI_Irecv(&receiveLeftBuffer, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &receiveRequest.back());
  }
  else
  {
    leftMessage = 1;
  }

  if (rank < numRanks - 1)
  {
    sendRequest.emplace_back();
    receiveRequest.emplace_back();

    MPI_Isend(&sendRightBuffer, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &sendRequest.back());
    MPI_Irecv(&receiveRightBuffer, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &receiveRequest.back());
  }
  else
  {
    rightMessage = 1;
  }

  MPI_Waitall(sendRequest.size(), sendRequest.data(), MPI_STATUS_IGNORE);
  MPI_Waitall(receiveRequest.size(), receiveRequest.data(), MPI_STATUS_IGNORE);

  if (rank > 0)
    leftMessage = receiveLeftBuffer;

  if (rank < numRanks - 1)
    rightMessage = receiveRightBuffer;

  const int sendReduceBuffer = leftMessage + rightMessage;
  int recvReduceBuffer;
  MPI_Allreduce(&sendReduceBuffer, &recvReduceBuffer, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0)
    std::cout << "Reduce is: " << recvReduceBuffer << " but is supposed to be: " << numRanks * 2 << std::endl;
}

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  // debugComm();
  // if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) == 0)
  // {
  //   int i=0;
  //   while (i == 0)
  //     sleep(5);
  // }

  ComputationParallel computationParallel = ComputationParallel();

  computationParallel.initialize(argc, argv);

  computationParallel.runSimulation();

  MPI_Finalize();
  return EXIT_SUCCESS;
}