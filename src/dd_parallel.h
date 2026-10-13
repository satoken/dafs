#ifndef DAFS_DD_PARALLEL_H
#define DAFS_DD_PARALLEL_H

#include <algorithm>
#include <array>
#include <exception>
#include <stdexcept>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace DDParallel {

// X folding, Y folding, and Z alignment own independent workspaces. As in
// IPknot's level parallelism, within-DP workers take precedence over this team.
inline int worker_limit(int threads, int dp_threads)
{
  if (threads < 0)
    throw std::invalid_argument("--dd-threads must be nonnegative");
#ifdef _OPENMP
  if (dp_threads != 1 || omp_in_parallel()) return 1;
  return std::min(3, threads == 0 ? omp_get_max_threads() : threads);
#else
  if (threads > 1)
    throw std::invalid_argument("--dd-threads > 1 requires OpenMP support");
  return 1;
#endif
}

// Join before aggregating objectives, updating multipliers or consuming any
// outputs. Exceptions must not escape an OpenMP region; rethrow in X/Y/Z order.
template <typename Decode>
int decode_subproblems(int threads, int dp_threads, Decode&& decode)
{
  const int limit = worker_limit(threads, dp_threads);
  if (limit == 1) {
    for (int job = 0; job < 3; ++job) decode(job);
    return 1;
  }
#ifdef _OPENMP
  std::array<std::exception_ptr, 3> errors;
  int workers = 1;
#pragma omp parallel num_threads(limit)
  {
#pragma omp single
    workers = omp_get_num_threads();
#pragma omp for schedule(static, 1)
    for (int job = 0; job < 3; ++job) {
      try { decode(job); }
      catch (...) { errors[job] = std::current_exception(); }
    }
  }
  for (const auto& error : errors)
    if (error) std::rethrow_exception(error);
  return workers;
#else
  return 1;
#endif
}

} // namespace DDParallel
#endif
