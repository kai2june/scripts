#pragma once
#include <vector>

template <typename T>
class ColumnMajorMatrix
{
  public:
    ColumnMajorMatrix(uint32_t dim1, uint32_t dim2) {}
  private:
    std::vector<std::vector<T>> all_column;
};