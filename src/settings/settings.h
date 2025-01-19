#pragma once

#include <iostream>
#include <array>

/** All settings that parametrize a simulation run.
 */
struct Settings
{
    std::array<int, 2> nCells;          //< number of cells in x and y direction
    std::array<double, 2> physicalSize; //< physical size of the domain
    double re = 1000;                   //< reynolds number
    double endTime = 10.0;              //< end time of the simulation
    double tau = 0.5;                   //< safety factor for time step width
    double maximumDt = 0.1;             //< maximum time step width

    std::array<double, 2> g{0., 0.}; //< external forces

    bool useDonorCell = false; //< if the donor cell scheme schould be used
    double alpha = 0.5;        //< factor for donor-cell scheme

    std::string pressureSolver = "SOR";  //< which pressure solver to use, "GaussSeidel" or "SOR"
    double omega = 1.0;                  //< overrelaxation factor
    double epsilon = 1e-5;               //< tolerance for the residual in the pressure solver
    int maximumNumberOfIterations = 1e5; //< maximum number of iterations in the solver

public:
    /// @brief Load and set parameters in struct from specified text file ! parse a text file with settings, each line contains "<parameterName> = <value>"
    /// @param filename Unix-type path to parameters file
    void loadParamsFromFile(std::string filename);

    //! output all settings to console
    void printSettings();

    /// @brief removes white space in the beginning of the line
    /// @param line line in file stream to be updated, will be mutated
    void removeWhitespaceAtBeginning(std::string &line);

private:

    /// @brief extracts the substring specifying the parameter name in the given file line
    /// @param line line in file stream to be operated on
    /// @return substring specifying the parameter name
    const std::string extractParameterName(std::string &line);

     /// @brief extract the substring specifying the parameter value in the given file line
    /// @param line line in file stream to be operated on
    /// @return substring specifying the parameter value
    const std::string extractValueString(std::string &line);
    
    /// @brief set specified parameter by autoconverting given string
    /// @param parameterName parameter identifier in struct to be changed
    /// @param valueString string of value that should be set, will be auto converted
    void setParameter(std::string &parameterName, std::string &valueString);

};