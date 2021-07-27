#pragma once
#include <vector>

template <typename T>
class RowMajorMatrix
{
  public:
    RowMajorMatrix(uint32_t dim1, uint32_t dim2) 
        : all_row(std::vector<std::vector<T>>(dim1, std::vector(dim2, 0))) 
        {}
  private:
    std::vector<std::vector<T>> all_row;
};