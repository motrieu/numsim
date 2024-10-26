# Master Lecture Numerical Simulation (University of Stuttgart)
Members: Moritz Trieu, Simon Grether, Nina Lock

## About
### Assignment 1
### Assignment 2
### Assignment 3
### Final Project

## How to Build
### Generating CMake
To generate all Cmakefiles, run either

``cmake`` in the root folder or alternatively run

``cmake ..`` in the /build/ folder

### Building Makefile
To build the entire project, run
``make`` in the /build/ folder

## How to use
The current entry point is ``numsim`` in /build/source/

## Testing
To run the test suit, run ``ctest`` in the /build/testing/ folder.

### Writing your own tests
Test can be found in the /testing/ folder. Here for every major module there is its own .cpp file. 
All test for that module can be found in this file. Each test file can therefore contain multiple suits.
If you create a new test file, make sure to include it in /testing/CMakeLists.txt as well as any needed dependecies,
mainly the modules tested in said test file.
