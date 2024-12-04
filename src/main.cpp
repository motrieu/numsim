#include "computation/computationParallel.h"

#include <array>
#include <mpi.h>

//#include <unistd.h>

int main(int argc, char *argv[])
{
  /*{
    int i=0;
    while (i == 0)
      sleep(5);
  }*/

  int rank;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_split(comm, rank < 48, rank, &comm);

  if (rank >= 48)
  {
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  ComputationParallel computationParallel = ComputationParallel();

  computationParallel.initialize(argc, argv);

  computationParallel.runSimulation();

  MPI_Finalize();

  return EXIT_SUCCESS;
}