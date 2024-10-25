#include "output_writer/write_paraview_output.h"
#include "settings/settings.h"

#include <iostream>
#include <cstdlib>

int main(int argc, char *argv[])
{
  // if the number of given command line arguments is only 1 (= the program name), print out usage information and exit
  if (argc == 1)
  {
    std::cout << "usage: " << argv[0] << " <filename>" << std::endl;
    return EXIT_FAILURE;
  }

  // read in the first argument
  std::string filename = argv[1];

  Settings settings;
  // load settings from file
  settings.loadFromFile(filename);
  // display all settings on console
  settings.printSettings();

  // write 5 output files
  for (int i = 0; i < 5; i++)
  {
    writeParaviewOutput(i);
  }

  return EXIT_SUCCESS;
}