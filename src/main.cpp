#include "computation/computationParallel.h"

#include <array>
#include <mpi.h>

//#include <unistd.h>

#include <chrono>
using namespace std::chrono;

int main(int argc, char *argv[])
{
  /*{
    int i=0;
    while (i == 0)
      sleep(5);
  }*/

  auto tSum = duration_cast<microseconds>(high_resolution_clock::now() - high_resolution_clock::now());

  std::cout << tSum.count() << std::endl;

  for (int i=0; i<5; i++){
    std::cout << i << std::endl;
    MPI_Init(&argc, &argv);

    ComputationParallel computationParallel = ComputationParallel();

    computationParallel.initialize(argc, argv);

    computationParallel.runSimulation();

    MPI_Finalize();

    tSum += duration_cast<microseconds>(high_resolution_clock::now() - tbegin);
  }

  std::cout << tSum.count() << std::endl;


  return EXIT_SUCCESS;
}