#include "computation/computation.h"
#include "partitioning/partitioning.h"

#include <array>

int main(int argc, char *argv[])
{

  std::array<int,2> nCellsGlobal{188,29};

  Partitioning partitioning = Partitioning();

  partitioning.initialize(nCellsGlobal);


  /*Computation computation = Computation();

  computation.initialize(argc, argv);

  computation.runSimulation();*/

  return EXIT_SUCCESS;
}