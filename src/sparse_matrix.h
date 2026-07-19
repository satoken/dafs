/*
 * Sparse row matrix used for Lagrange multipliers.
 */

#ifndef __INC_SPARSE_MATRIX_H__
#define __INC_SPARSE_MATRIX_H__

#include "typedefs.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <unordered_map>
#include <vector>

class SparseFloatMatrix
{
public:
  SparseFloatMatrix() : columns_(0), nonzeros_(0) { }

  void assign(uint rows, uint columns)
  {
    rows_.assign(rows, Row());
    ordered_rows_.assign(rows, SV());
    ordered_rows_dirty_.assign(rows, false);
    columns_ = columns;
    nonzeros_ = 0;
  }

  uint rows() const { return rows_.size(); }
  uint columns() const { return columns_; }
  size_t nonzeros() const { return nonzeros_; }

  void add(uint row, uint column, float value)
  {
    if (value == 0.0f)
      return;
    set(row, column, get(row, column) + value);
  }

  void scale(float factor)
  {
    if (factor == 0.0f) {
      for (uint row = 0; row < rows_.size(); ++row) {
        rows_[row].clear();
        ordered_rows_[row].clear();
        ordered_rows_dirty_[row] = false;
      }
      nonzeros_ = 0;
      return;
    }
    for (uint row = 0; row < rows_.size(); ++row) {
      for (auto& entry : rows_[row])
        entry.second *= factor;
      ordered_rows_dirty_[row] = true;
    }
  }

  // Remove small entries and optionally cap probabilities.  Iterating only
  // over materialized entries keeps thresholding proportional to nnz.
  void prune(float cutoff, float maximum = 0.0f)
  {
    for (uint row = 0; row < rows_.size(); ++row) {
      Row& entries = rows_[row];
      for (auto it = entries.begin(); it != entries.end(); ) {
        if (it->second <= cutoff) {
          it = entries.erase(it);
          --nonzeros_;
        } else {
          if (maximum > 0.0f && it->second > maximum)
            it->second = maximum;
          ++it;
        }
      }
      ordered_rows_dirty_[row] = true;
    }
  }

  VVF dense() const
  {
    VVF matrix(rows(), VF(columns_, 0.0f));
    for (uint row = 0; row < rows_.size(); ++row)
      for (const auto& [column, value] : rows_[row])
        matrix[row][column] = value;
    return matrix;
  }

  // Sorted row view for decoders that visit columns monotonically.  It is
  // rebuilt only after that row has been updated.
  const SV& ordered_row(uint row) const
  {
    assert(row < rows_.size());
    if (ordered_rows_dirty_[row]) {
      SV& ordered = ordered_rows_[row];
      ordered.assign(rows_[row].begin(), rows_[row].end());
      std::sort(ordered.begin(), ordered.end(),
                [](const auto& lhs, const auto& rhs) {
                  return lhs.first < rhs.first;
                });
      ordered_rows_dirty_[row] = false;
    }
    return ordered_rows_[row];
  }

  float get(uint row, uint column) const
  {
    assert(row < rows_.size() && column < columns_);
    const auto it = rows_[row].find(column);
    return it == rows_[row].end() ? 0.0f : it->second;
  }

  void set(uint row, uint column, float value)
  {
    assert(row < rows_.size() && column < columns_);
    Row& entries = rows_[row];
    const auto it = entries.find(column);
    if (value == 0.0f) {
      if (it != entries.end()) {
        entries.erase(it);
        --nonzeros_;
        ordered_rows_dirty_[row] = true;
      }
    } else if (it == entries.end()) {
      entries.emplace(column, value);
      ++nonzeros_;
      ordered_rows_dirty_[row] = true;
    } else {
      it->second = value;
      ordered_rows_dirty_[row] = true;
    }
  }

private:
  using Row = std::unordered_map<uint, float>;
  std::vector<Row> rows_;
  mutable std::vector<SV> ordered_rows_;
  mutable std::vector<bool> ordered_rows_dirty_;
  uint columns_;
  size_t nonzeros_;
};

#endif // __INC_SPARSE_MATRIX_H__
