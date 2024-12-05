#pragma once

#include "storage/fieldVariable.h"
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

    const FieldVariable& pOld() const;
    const FieldVariable& pOld2() const;
 	

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

    /// @brief get constant value of p in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return constant value of p in element with indices i,j
    double pOld(int i, int j) const;
 
    /// @brief get reference to value of p in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of p in element with indices i,j
    double&	pOld(int i, int j);

    /// @brief get constant value of p in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return constant value of p in element with indices i,j
    double pOld2(int i, int j) const;
 
    /// @brief get reference to value of p in element i,j
    /// @param i index of element in x-direction
    /// @param j index of element in y-direction
    /// @return reference to value of p in element with indices i,j
    double&	pOld2(int i, int j);

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

    FieldVariable pOld_;
    FieldVariable pOld2_;

    /// @brief stores the values of the preliminary velocity in x-direction, lives on the right face of each cell
    FieldVariable f_;

    /// @brief stores the values of the preliminary velocity in y-direction, lives on the upper face of each cell
    FieldVariable g_;

    /// @brief stores the values of the right hand side, lives in the centre of each cell
    FieldVariable rhs_;

};