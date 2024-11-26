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

  MPI_Init(&argc, &argv);

  ComputationParallel computationParallel = ComputationParallel();

  computationParallel.initialize(argc, argv);

  computationParallel.runSimulation();

  MPI_Finalize();

  return EXIT_SUCCESS;
}