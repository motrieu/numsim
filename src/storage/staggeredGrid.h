#pragma once

#include "fieldVariable.h"
#include <array>

class StaggeredGrid
{
public:
   
    /// @brief constructor of staggered grid
    /// @param nCells two-dimensional array for number of elements in x and y direction
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
 

    /// @brief get mesh width in x-direction
    /// @return mesh width in x-direction
    double dx() const;
    
    /// @brief get mesh width in y-direction
    /// @return mesh width in y-direction
    double dy() const;
 
    
    /// @brief get first valid index for u in x-direction
    /// @return first valid index for u in x-direction
    int	uIBegin() const;
    
    /// @brief get one after last valid index for u in x-direction
    /// @return one after last valid index for u in x-direction
    int	uIEnd() const;
    
    /// @brief get first valid index for u in y-direction
    /// @return first valid index for u in y-direction
    int	uJBegin() const;
 
    /// @brief get one after last valid index for u in y-direction
    /// @return one after last valid index for u in y-direction
    int uJEnd() const;
    
    /// @brief get first valid index for v in x-direction
    /// @return first valid index for v in x-direction
    int	vIBegin() const;
    
    /// @brief get one after last valid index for v in x-direction
    /// @return one after last valid index for v in x-direction
    int	vIEnd() const;
    
    /// @brief get first valid index for v in y-direction
    /// @return first valid index for v in y-direction
    int	vJBegin() const;

    /// @brief get one after last valid index for v in y-direction
    /// @return one after last valid index for v in y-direction
    int	vJEnd() const;
    
    /// @brief get first valid index for p in x-direction
    /// @return first valid index for p in x-direction
    int	pIBegin() const;
    
    /// @brief get one after last valid index for p in x-direction
    /// @return one after last valid index for p in x-direction
    int	pIEnd() const;

    /// @brief get first valid index for p in y-direction
    /// @return first valid index for p in y-direction
    int	pJBegin() const;
    
    /// @brief get one after last valid index for p in y-direction
    /// @return one after last valid index for p in y-direction
    int	pJEnd() const;

protected:
    const std::array<int,2> nCells_;
    const std::array<double,2> meshWidth_;
 
    FieldVariable u_;
    FieldVariable v_;
    FieldVariable p_;
    FieldVariable rhs_;
    FieldVariable f_;
    FieldVariable g_;

};