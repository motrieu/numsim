#include "array2dInt.h"

#include <cassert>

Array2DInt::Array2DInt(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to -1
  data_.resize(size_[0]*size_[1], -1);
}

//! get the size
std::array<int,2> Array2DInt::size() const
{
  return size_;
}

int &Array2DInt::operator()(int i, int j)
{
  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}

int Array2DInt::operator()(int i, int j) const
{
  const int index = j*size_[0] + i;

  // assert that indices are in range
  assert(0 <= i && i < size_[0]);
  assert(0 <= j && j < size_[1]);
  assert(j*size_[0] + i < (int)data_.size());

  return data_[index];
}