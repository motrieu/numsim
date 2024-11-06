#include "computation/computation.h"

int main(int argc, char *argv[])
{

  Computation computation = Computation();

  computation.initialize(argc, argv);

  computation.runSimulation();

  return EXIT_SUCCESS;
}