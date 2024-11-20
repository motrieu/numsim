#include "computation/computation.h"

#include <chrono>
using namespace std::chrono;

int main(int argc, char *argv[])
{

  auto tSum = duration_cast<microseconds>(high_resolution_clock::now() - high_resolution_clock::now());

  std::cout << tSum.count() << std::endl;

  for (int i=0; i<1; i++){
    std::cout << i << std::endl;
    auto tbegin = high_resolution_clock::now();

    Computation computation = Computation();

    computation.initialize(argc, argv);

    computation.runSimulation();

    tSum += duration_cast<microseconds>(high_resolution_clock::now() - tbegin);
  }

  std::cout << tSum.count() << std::endl;
  

  return EXIT_SUCCESS;
}