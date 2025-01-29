#pragma once

#include "storage/fieldVariable.h"
#include "storage/array2d.h"
#include "storage/array2dInt.h"
#include <array>

class StaggeredGrid
{
public:
   
    /// @brief constructor of staggered grid
    /// @param nCells two-dimensional array for number of elements in x and y direction (halo cells not included)
    /// @param meshWidth two-dimensional array for mesh width in x and y direction
    StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth);

    /// @brief returns mesh width of staggered grid in x and y direction
    /// @return two-dimensional array for mesh width in x and y direction
    const std::array<double,2> meshWidth() const;

    /// @brief returns number of elements in x and y direction
    /// @return two-dimensional array for number of elements in x and y direction
    const std::array<int,2> nCells() const;

    /// @brief  get field variable u
    /// @return reference to field variable u
    const FieldVariable& u() const;

    /// @brief  get field variable v
    /// @return reference to field variable v
    const FieldVariable& v() const;

    /// @brief  get field variable p
    /// @return reference to field variable p
    const FieldVariable& p() const;
 	

    /// @brief get constant value of u in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return constant value of u in element with indices i,j
    double u(int i, int j) const;
 
    /// @brief get reference to value of u in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of u in element with indices i,j
    double& u(int i, int j);
 
    /// @brief get constant value of v in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return constant value of v in element with indices i,j
    double v(int i, int j) const;
 
    /// @brief get reference to value of v in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of v in element with indices i,j
    double& v(int i, int j);
 
    /// @brief get constant value of p in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return constant value of p in element with indices i,j
    double p(int i, int j) const;
 
    /// @brief get reference to value of p in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of p in element with indices i,j
    double&	p(int i, int j);

    /// @brief get reference to value of the rhs in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of the rhs in element with indices i,j
    double&	rhs(int i, int j);
 
    /// @brief get reference to value of F in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of F in element with indices i,j
    double&	f(int i, int j);
 
    /// @brief get reference to value of G in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of G in element with indices i,j
    double& g(int i, int j);

    /// @brief get constant setup value of element i,j 
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return setup value of element i,j
    int setup(int i, int j) const;

    /// @brief get reference to setup value of element i,j 
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to setup value of element i,j
    int& setup(int i, int j);

    /// @brief 
    /// @param x 
    /// @param y 
    /// @return 
    int interpolationContainsNoSlip(double x, double y) const;

    /// @brief get constant edgeDirection value of element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return edgeDirection value of element i,j
    int edgeDirections(int i, int j) const;
    
    /// @brief get reference to edgeDirection value of element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to edgeDirection value of element i,j
    int& edgeDirections(int i, int j);


    /// @brief get constant numberFaces value of element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return numberFaces value of element i,j
    int numberFaces(int i, int j) const;

    /// @brief get reference to numberFaces value of element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to numberFaces value of element i,j
    int& numberFaces(int i, int j);

    /// @brief get constant Dirichlet boundary value uIn for element i,j used for INFLOW condition
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @return Dirichlet boundary value uIn for element i,j
    double uIn(int i, int j) const;

    /// @brief get reference to Dirichlet boundary value uIn for element i,j used for INFLOW condition
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @return reference to Dirichlet boundary value uIn for element i,j
    double& uIn(int i, int j);

    /// @brief get constant Dirichlet boundary value vIn for element i,j used for INFLOW condition
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @return Dirichlet boundary value vIn for element i,j
    double vIn(int i, int j) const;

    /// @brief get reference to Dirichlet boundary value vIn for element i,j used for INFLOW condition
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @return reference to Dirichlet boundary value vIn for element i,j
    double& vIn(int i, int j);

    /// @brief get constant Dirichlet boundary value pRB for element i,j used for PRESSURE condition
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @return Dirichlet boundary value vIn for element i,j
    double pRB(int i, int j) const;

    /// @brief get reference to Dirichlet boundary value pRB for element i,j used for PRESSURE condition
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction 
    /// @return reference to Dirichlet boundary value vIn for element i,j
    double& pRB(int i, int j);
 
    
    /// @brief index that indicates a fluid cell
    /// @return 0 
    int indexFluid() const;

    /// @brief index that indicates an obstacle cell with NOSLIP condition
    /// @return 1 
    int indexNoSlip() const;

    /// @brief index that indicates an obstacle cell with SLIP condition
    /// @return 2
    int indexSlip() const;

    /// @brief index that indicates an obstacle cell with INFLOW condition, combined with Pressure Neumann
    /// @return  
    int indexInflow() const;

    /// @brief index that indicates an obstacle cell with OUTFLOW condition, combined with Pressure Neumann
    /// @return 4
    int indexOutflow() const;

    /// @brief index that indicates an obstacle cell with PRESSURE (Pressure Dirichlet) condition, combined with OUTFLOW 
    /// @return 5
    int indexPressure() const;

    /// @brief get mesh width in x-direction
    /// @return mesh width in x-direction
    double dx() const;
    
    /// @brief get mesh width in y-direction
    /// @return mesh width in y-direction
    double dy() const;
 
    
    /// @brief get first inner index for u in x-direction
    /// @return first inner index for u in x-direction
    int	uIBegin() const;
    
    /// @brief get one after last inner index for u in x-direction
    /// @return one after last inner index for u in x-direction
    int	uIEnd() const;
    
    /// @brief get first inner index for u in y-direction
    /// @return first inner index for u in y-direction
    int	uJBegin() const;
 
    /// @brief get one after last inner index for u in y-direction
    /// @return one after last inner index for u in y-direction
    int uJEnd() const;
    
    /// @brief get first inner index for v in x-direction
    /// @return first inner index for v in x-direction
    int	vIBegin() const;
    
    /// @brief get one after last inner index for v in x-direction
    /// @return one after last inner index for v in x-direction
    int	vIEnd() const;
    
    /// @brief get first inner index for v in y-direction
    /// @return first inner index for v in y-direction
    int	vJBegin() const;

    /// @brief get one after last inner index for v in y-direction
    /// @return one after last inner index for v in y-direction
    int	vJEnd() const;
    
    /// @brief get first inner index for p in x-direction
    /// @return first inner index for p in x-direction
    int	pIBegin() const;
    
    /// @brief get one after last inner index for p in x-direction
    /// @return one after last inner index for p in x-direction
    int	pIEnd() const;

    /// @brief get first inner index for p in y-direction
    /// @return first inner index for p in y-direction
    int	pJBegin() const;
    
    /// @brief get one after last inner index for p in y-direction
    /// @return one after last inner index for p in y-direction
    int	pJEnd() const;

    
    /// @brief get first valid index for setup in x-direction
    /// @return 0
    int	setupIBegin() const;

    /// @brief get last valid index for setup in x-direction
    /// @return nCells_[0] + 2
    int	setupIEnd() const;
    
    /// @brief get first valid index for setup in y-direction
    /// @return 0
    int	setupJBegin() const;
    
    /// @brief get last valid index for setup in y-direction
    /// @return nCells_[1] + 2
    int setupJEnd() const;

protected:
    /// @brief two-dimensional array for number of elements in x and y direction (halo cells not included)
    const std::array<int,2> nCells_;

    /// @brief two-dimensional array for mesh width in x and y direction
    const std::array<double,2> meshWidth_;
 
    /// @brief stores the values of the velocity in x-direction, lives on the right face of each cell
    FieldVariable u_;

    /// @brief stores the values of the velocity in y-direction, lives on the upper face of each cell
    FieldVariable v_;

    /// @brief stores the values of the pressure, lives in the centre of each cell
    FieldVariable p_;

    /// @brief stores the values of the preliminary velocity in x-direction, lives on the right face of each cell
    FieldVariable f_;

    /// @brief stores the values of the preliminary velocity in y-direction, lives on the upper face of each cell
    FieldVariable g_;

    /// @brief stores the values of the right hand side, lives in the centre of each cell
    FieldVariable rhs_;

    /// @brief stores the setup values for all elements
    ///        setup describes whether a cell is a fluid cell or an obstacle cell
    ///        furthermore, if it is an obstacle cell, setup determines which boundary condition holds for this obstacle
    ///        0 for fluid cell, 1 for obstacle cell and NOSLIP, 2 for obstacle cell and SLIP, 3 for obstacle cell and INFLOW, 4 for obstacle cell and OUTFLOW, 5 for obstacle cell and PRESSURE
    ///        default: -1 (gets overwritten in every cell)
    Array2DInt setup_;

    /// @brief stores the edgeDirections for all elements
    ///        edgeDirections describe the direction in which the obstacle element has a fluid neighbor (default: -1 if it has none)
    ///        0 for north, 1 for north-east corner, 2 for east, 3 for east-south corner, 4 for south, 5 for south-west corner, 6 for west, 7 for west-north corner
    Array2DInt edgeDirections_;

    /// @brief stores the numberFaces for all elements
    ///        for corner fluid cells: numberFaces = 2
    ///        for diagonal fluid cells: numberFaces = 0
    ///        numberFaces is needed to distinguish between corner fluid cells and diagonal fluid cells
    ///        corner fluid cells: obstacle has two neighbouring fluid cells building a "fluid corner"
    ///        diagonal fluid cells: obstacle has no neighbouring fluid cell, but its diagonal cell is a fluid cell
    ///        numberFaces is default -1 (gets overwritten in every cell)
    Array2DInt numberFaces_;

    /// @brief stores the Dirichlet boundary values of uIn for all elements (needed for INFLOW condition)
    Array2D uIn_;

    /// @brief stores the Dirichlet boundary values of vIn for all elements (needed for INFLOW condition)
    Array2D vIn_;

    /// @brief stores the Dirichlet boundary values of pRB for all elements (needed for PRESSURE condition)
    Array2D pRB_;

};