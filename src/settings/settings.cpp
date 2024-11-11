#include "settings.h"
#include <fstream>
#include <iomanip>

void Settings::loadFromFile(std::string filename)
{
    // open file
    std::ifstream file(filename.c_str(), std::ios::in);

    // check if file is open
    if (!file.is_open())
    {
        std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
        return;
    }

    // loop over lines of file
    for (int lineNo = 0;; lineNo++)
    {
        // read line
        std::string line;
        getline(file, line);

        // at the end of the file break for loop
        if (file.eof())
            break;

        removeWhitespaceAtBeginning(line);

        std::string parameterName;
        std::string valueString;

        // checks if line is a comment and if line contains a "="-sign
        // if it contains a "="-sign (i.e. a parameter is defined) and it is not a comment, the parameter is parsed, otherwise the line is ignored
        if ((line[0] != '#') && (line.find("=") != std::string::npos))
        {
            parameterName = extractParameterName(line);
            valueString = extractValueString(line);

            setParameter(parameterName, valueString);
        }
    }
}

void Settings::setParameter(std::string &parameterName, std::string &valueString)
{
    if (parameterName == "physicalSizeX")
        physicalSize[0] = std::stod(valueString);
    else if (parameterName == "physicalSizeY")
        physicalSize[1] = std::stod(valueString);
    else if (parameterName == "endTime")
        endTime = std::stod(valueString);
    else if (parameterName == "re")
        re = std::stod(valueString);
    else if (parameterName == "gX")
        g[0] = std::stod(valueString);
    else if (parameterName == "gY")
        g[1] = std::stod(valueString);
    else if (parameterName == "dirichletBottomX")
        dirichletBcBottom[0] = std::stod(valueString);
    else if (parameterName == "dirichletBottomY")
        dirichletBcBottom[1] = std::stod(valueString);
    else if (parameterName == "dirichletTopX")
        dirichletBcTop[0] = std::stod(valueString);
    else if (parameterName == "dirichletTopY")
        dirichletBcTop[1] = std::stod(valueString);
    else if (parameterName == "dirichletLeftX")
        dirichletBcLeft[0] = std::stod(valueString);
    else if (parameterName == "dirichletLeftY")
        dirichletBcLeft[1] = std::stod(valueString);
    else if (parameterName == "dirichletRightX")
        dirichletBcRight[0] = std::stod(valueString);
    else if (parameterName == "dirichletRightY")
        dirichletBcRight[1] = std::stod(valueString);
    else if (parameterName == "nCellsX")
        nCells[0] = (int)std::stod(valueString) + 2;
    else if (parameterName == "nCellsY")
        nCells[1] = (int)std::stod(valueString) + 2;
    else if (parameterName == "useDonorCell")
    {
        if (valueString == "true")
            useDonorCell = true;
        else if (valueString == "false")
            useDonorCell = false;
        else
            throw std::invalid_argument("Assigned value for useDonorCell is not of type boolean.");
    }
    else if (parameterName == "alpha")
        alpha = std::stod(valueString);
    else if (parameterName == "tau")
        tau = std::stod(valueString);
    else if (parameterName == "maximumDt")
        maximumDt = std::stod(valueString);
    else if (parameterName == "pressureSolver")
        pressureSolver = valueString;
    else if (parameterName == "omega")
        omega = std::stod(valueString);
    else if (parameterName == "epsilon")
        epsilon = std::stod(valueString);
    else if (parameterName == "maximumNumberOfIterations")
        maximumNumberOfIterations = (int)std::stod(valueString);
    else
        throw std::invalid_argument("The parameter " + parameterName + " is not implemented.");
}

const std::string Settings::extractValueString(std::string &line)
{
    std::string extractedString;
    int valueStartIndex = line.find_first_not_of(" \t", line.find("=") + 1);
    extractedString = line.substr(valueStartIndex, line.find_first_of(" #\t\n", valueStartIndex) - valueStartIndex);

    return extractedString;
}

const std::string Settings::extractParameterName(std::string &line)
{
    std::string paramName;
    paramName = line.substr(0, line.find_first_of(" =\t"));
    return paramName;
}
void Settings::removeWhitespaceAtBeginning(std::string &line)
{
    line.erase(0, line.find_first_not_of(" \t"));
}

void Settings::printSettings()
{
    std::cout << "Settings: " << std::endl
              << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
              << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
              << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"
              << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"
              << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
              << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
              << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
              << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}