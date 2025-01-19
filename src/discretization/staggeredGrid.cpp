#include "staggeredGrid.h"

StaggeredGrid::StaggeredGrid(std::array<int,2> nCells, std::array<double,2> meshWidth) :
    nCells_(nCells), 
    meshWidth_(meshWidth), 
    u_(FieldVariable({nCells[0]+2, nCells[1]+2}, {0.0, -0.5*meshWidth[1]}, meshWidth)), // first (i,j=0) u node (halo node) lives half a y-mesh width below the cartesian origin (x,y=0)
    v_(FieldVariable({nCells[0]+2, nCells[1]+2}, {-0.5*meshWidth[0], 0.0}, meshWidth)), // first (i,j=0) v node (halo node) lives half a x-mesh width left to the cartesian origin (x,y=0)
    p_(FieldVariable({nCells[0]+2, nCells[1]+2}, {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)), // first (i,j=0) p node (halo node) lives half a x- and half a y-mesh width apart from the cartesian origin (x,y=0)
    f_(FieldVariable({nCells[0]+2, nCells[1]+2}, {0.0, -0.5*meshWidth[1]}, meshWidth)), // first (i,j=0) f node (halo node) lives half a y-mesh width below the cartesian origin (x,y=0) (analog to u)
    g_(FieldVariable({nCells[0]+2, nCells[1]+2}, {-0.5*meshWidth[0], 0.0}, meshWidth)), // first (i,j=0) g node (halo node) lives half a x-mesh width left to the cartesian origin (x,y=0) (analog to v)
    rhs_(FieldVariable({nCells[0]+2, nCells[1]+2}, {-0.5*meshWidth[0], -0.5*meshWidth[1]}, meshWidth)), // first (i,j=0) rhs node (halo node) lives half a x- and half a y-mesh width apart from the cartesian origin (x,y=0) (analog to p)
    setup_({nCells[0]+2, nCells[1]+2}),
    edgeDirections_({nCells[0]+2, nCells[1]+2}),
    uIn_({nCells[0]+2, nCells[1]+2}),
    vIn_({nCells[0]+2, nCells[1]+2}),
    pRB_({nCells[0]+2, nCells[1]+2})
{
}

const std::array<double,2> StaggeredGrid::meshWidth() const
{
    return meshWidth_;
}

const std::array<int,2> StaggeredGrid::nCells() const
{
    return nCells_;
}

const FieldVariable& StaggeredGrid::u() const
{
    return u_;
}

const FieldVariable& StaggeredGrid::v() const
{
    return v_;
}

const FieldVariable& StaggeredGrid::p() const
{
    return p_;
}

double StaggeredGrid::u(int i, int j) const
{
    return u_(i,j);
}
 
double& StaggeredGrid::u(int i, int j)
{
    return u_(i,j);
}
 
double StaggeredGrid::v(int i, int j) const
{
    return v_(i,j);
}
 
double& StaggeredGrid::v(int i, int j)
{
    return v_(i,j);
}
 
double StaggeredGrid::p(int i, int j) const
{
    return p_(i,j);
}
 
double& StaggeredGrid::p(int i, int j)
{
    return p_(i,j);
}

double&	StaggeredGrid::rhs(int i, int j)
{
    return rhs_(i,j);
}
 
double&	StaggeredGrid::f(int i, int j)
{
    return f_(i,j);
}
 
double& StaggeredGrid::g(int i, int j)
{
    return g_(i,j);
}


int StaggeredGrid::setup(int i, int j) const
{
    return setup_(i,j);
}

int& StaggeredGrid::setup(int i, int j)
{
    return setup_(i,j);
}

int StaggeredGrid::edgeDirections(int i, int j) const
{
    return edgeDirections_(i,j);
}

int& StaggeredGrid::edgeDirections(int i, int j)
{
    return edgeDirections_(i,j);
}

double StaggeredGrid::uIn(int i, int j) const
{
    return uIn_(i,j);
}

double& StaggeredGrid::uIn(int i, int j)
{
    return uIn_(i,j);
}

double StaggeredGrid::vIn(int i, int j) const
{
    return vIn_(i,j);
}

double& StaggeredGrid::vIn(int i, int j)
{
    return vIn_(i,j);
}

double StaggeredGrid::pRB(int i, int j) const
{
    return pRB_(i,j);
}

double& StaggeredGrid::pRB(int i, int j)
{
    return pRB_(i,j);
}

int StaggeredGrid::indexFluid()
{
    return 0;
}

int StaggeredGrid::indexNoSlip()
{
    return 1;
}

int StaggeredGrid::indexSlip()
{
    return 2;
}

int StaggeredGrid::indexInflow()
{
    return 3;
}

int StaggeredGrid::indexOutflow()
{
    return 4;
}

double StaggeredGrid::dx() const 
{
    return meshWidth_[0];
}
 
double StaggeredGrid::dy() const
{
    return meshWidth_[1];
}

int	StaggeredGrid::uIBegin() const 
{
    return 1;
}
 
int	StaggeredGrid::uIEnd() const
{
    return nCells_[0];
}
 
int	StaggeredGrid::uJBegin() const
{
    return 1;
}
 
int StaggeredGrid::uJEnd() const
{
    return nCells_[1] + 1;
}

int	StaggeredGrid::vIBegin() const 
{
    return 1;
}
 
int	StaggeredGrid::vIEnd() const
{
    return nCells_[0] + 1;
}
 
int	StaggeredGrid::vJBegin() const
{
    return 1;
}
 
int StaggeredGrid::vJEnd() const
{
    return nCells_[1];
}

int	StaggeredGrid::pIBegin() const 
{
    return 1;
}
 
int	StaggeredGrid::pIEnd() const
{
    return nCells_[0] + 1;
}
 
int	StaggeredGrid::pJBegin() const
{
    return 1;
}
 
int StaggeredGrid::pJEnd() const
{
    return nCells_[1] + 1;
}

int	StaggeredGrid::setupIBegin() const 
{
    return 0;
}
 
int	StaggeredGrid::setupIEnd() const
{
    return nCells_[0] + 2;
}
 
int	StaggeredGrid::setupJBegin() const
{
    return 0;
}
 
int StaggeredGrid::setupJEnd() const
{
    return nCells_[1] + 2;
}
