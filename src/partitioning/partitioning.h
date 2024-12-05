#pragma once

#include <array>
#include <vector>
#include <cmath>
#include <cassert>
#include <mpi.h>

class Partitioning
{
public:
  //! compute partitioning, set internal variables
  //! calculates how the domain is decomposed into partitions based on number of ranks
  void initialize(std::array<int, 2> nCellsGlobal);

  //! get the local number of cells in the own subdomain
  std::array<int, 2> nCellsLocal() const;

  //! get the global number of cells in the whole computational domain
  //! used in OutputWriterParaviewParallel
  std::array<int, 2> nCellsGlobal() const;

  //! get the own MPI rank no
  //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
  int ownRankNo() const;

  //! number of MPI ranks
  int nRanks() const;

  //! if the own partition has part of the bottom boundary of the whole domain
  bool ownPartitionContainsBottomBoundary() const;

  //! if the own partition has part of the top boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsTopBoundary() const;

  //! if the own partition has part of the left boundary of the whole domain
  bool ownPartitionContainsLeftBoundary() const;

  //! if the own partition has part of the right boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsRightBoundary() const;

  //! get the rank no of the left neighbouring rank
  int leftNeighbourRankNo() const;

  //! get the rank no of the right neighbouring rank
  int rightNeighbourRankNo() const;

  //! get the rank no of the top neighbouring rank
  int topNeighbourRankNo() const;

  //! get the rank no of the bottom neighbouring rank
  int bottomNeighbourRankNo() const;

  //! get the offset values for counting local nodes in x and y direction.
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  //! used in OutputWriterParaviewParallel
  std::array<int, 2> nodeOffset() const;

private:

  /// @brief rounds double value to nearest integer
  /// @param x double value
  /// @return rounded integer value
  int roundToInt(double x);


  /// @brief number of partitions/ranks, defined by how mpi has been executed
  int numberOfRanks_;


  /// @brief index of the mpi rank of this process
  int ownRank_;


  /// @brief index of the mpi rank of the left neighbour of this process,
  ///        is -1 if this process has no left neighbour, i.e. has a left Dirichlet boundary
  int rankLeft_;
  
  /// @brief index of the mpi rank of the right neighbour of this process,
  ///        is -1 if this process has no right neighbour, i.e. has a right Dirichlet boundary
  int rankRight_;
  
  /// @brief index of the mpi rank of the lower neighbour of this process,
  ///        is -1 if this process has no lower neighbour, i.e. has a bottom Dirichlet boundary
  int rankLower_;
  
  /// @brief index of the mpi rank of the upper neighbour of this process,
  ///        is -1 if this process has no upper neighbour, i.e. has a top Dirichlet boundary
  int rankUpper_;


  /// @brief is true if this process has a Dirichlet boundary on the left
  bool ownPartitionContainsLeftBoundary_;

  /// @brief is true if this process has a Dirichlet boundary on the right
  bool ownPartitionContainsRightBoundary_;

  /// @brief is true if this process has a Dirichlet boundary at the bottom
  bool ownPartitionContainsBottomBoundary_;

  /// @brief is true if this process has a Dirichlet boundary at the top
  bool ownPartitionContainsTopBoundary_;

  /// @brief two-dimensional array for number of elements in x and y direction for this process (halo cells not included)
  std::array<int,2> nCellsLocal_;

  /// @brief two-dimensional array for number of elements in x and y direction for global physical domain (halo cells not included)
  std::array<int,2> nCellsGlobal_;

  /// @brief offset values for counting local nodes in x and y direction
  ///        (i_local,j_local) + nodeOffset = (i_global,j_global)
  ///        sed in OutputWriterParaviewParallel
  std::array<int,2> nodeOffset_;

};
