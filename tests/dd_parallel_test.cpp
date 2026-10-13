#include "dd_parallel.h"

#include <array>
#include <iostream>
#include <set>
#include <string>

namespace {
void require(bool value, const char* message)
{
  if (!value) throw std::runtime_error(message);
}

void assignments(int requested, int dp_threads, int expected)
{
  std::array<int, 3> calls{}, owners{}, values{};
  const int workers = DDParallel::decode_subproblems(requested, dp_threads,
      [&](int job) {
    ++calls[job];
#ifdef _OPENMP
    owners[job] = omp_get_thread_num();
    require(omp_get_level() <= 1, "nested DD worker team");
#endif
    values[job] = 10 + job;
  });
  require(workers == expected, "wrong number of DD workers");
  require(std::set<int>(owners.begin(), owners.end()).size() == size_t(expected),
          "subproblems were not assigned to distinct workers");
  require(calls == std::array<int, 3>{1, 1, 1}, "missing or repeated subproblem");
  require(values == std::array<int, 3>{10, 11, 12}, "outputs used before joining");
}

void exceptions(int threads)
{
  std::array<int, 3> calls{};
  bool caught = false;
  try {
    DDParallel::decode_subproblems(threads, 1, [&](int job) {
      ++calls[job];
      if (job != 2) throw std::runtime_error(job == 0 ? "X failure" : "Y failure");
    });
  } catch (const std::runtime_error& error) {
    caught = std::string(error.what()) == "X failure";
  }
  require(caught, "DD exception lost or rethrown in nondeterministic order");
  require(calls == (threads == 1 ? std::array<int, 3>{1, 0, 0}
                                : std::array<int, 3>{1, 1, 1}),
          "DD worker exceptions were not joined");
}
}

int main()
{
  try {
    bool rejected = false;
    try { DDParallel::worker_limit(-1, 1); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "negative DD threads accepted");
    assignments(1, 1, 1);
    exceptions(1);
#ifdef _OPENMP
    const int old_dynamic = omp_get_dynamic();
    const int old_threads = omp_get_max_threads();
    omp_set_dynamic(0);
    omp_set_num_threads(3);
    assignments(0, 1, 3);
    assignments(2, 1, 2);
    assignments(3, 1, 3);
    assignments(8, 1, 3);
    assignments(3, 0, 1);
    assignments(3, 2, 1);
    exceptions(3);
    int nested_workers = 0;
#pragma omp parallel num_threads(2)
    {
#pragma omp single
      nested_workers = DDParallel::decode_subproblems(3, 1, [](int) {
        require(omp_get_level() == 1, "nested worker team created");
      });
    }
    require(nested_workers == 1, "nested DD was not serialized");
    omp_set_num_threads(old_threads);
    omp_set_dynamic(old_dynamic);
#else
    assignments(0, 1, 1);
    rejected = false;
    try { DDParallel::worker_limit(2, 1); }
    catch (const std::invalid_argument&) { rejected = true; }
    require(rejected, "DD parallelism accepted without OpenMP");
#endif
    std::cout << "DD worker ownership, join, nesting and exceptions passed\n";
  } catch (const std::exception& error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
