#include "computation/computationParallel.h"

#include <array>

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);

  ComputationParallel computationParallel = ComputationParallel();

  computationParallel.initialize(argc, argv);

  computationParallel.runSimulation();

  MPI_Finalize();

  return EXIT_SUCCESS;
}