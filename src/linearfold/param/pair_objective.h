#pragma once

#include <cstddef>
#include <string>

// Zero-energy parameter model for maximizing an externally supplied
// base-pair objective with the LinearFold Viterbi grammar.
class PairObjectiveNearestNeighbor
{
public:
    using ScoreType = float;

    explicit PairObjectiveNearestNeighbor(const std::string&) { }

    ScoreType score_hairpin(size_t, size_t) const { return 0.0f; }
    ScoreType score_single_loop(size_t, size_t, size_t, size_t) const { return 0.0f; }
    ScoreType score_helix(size_t, size_t, size_t) const { return 0.0f; }
    ScoreType score_multi_loop(size_t, size_t) const { return 0.0f; }
    ScoreType score_multi_paired(size_t, size_t) const { return 0.0f; }
    ScoreType score_multi_unpaired(size_t, size_t) const { return 0.0f; }
    ScoreType score_external_zero() const { return 0.0f; }
    ScoreType score_external_paired(size_t, size_t) const { return 0.0f; }
    ScoreType score_external_unpaired(size_t, size_t) const { return 0.0f; }

    void count_hairpin(size_t, size_t, ScoreType) { }
    void count_single_loop(size_t, size_t, size_t, size_t, ScoreType) { }
    void count_helix(size_t, size_t, size_t, ScoreType) { }
    void count_multi_loop(size_t, size_t, ScoreType) { }
    void count_multi_paired(size_t, size_t, ScoreType) { }
    void count_multi_unpaired(size_t, size_t, ScoreType) { }
    void count_external_zero(ScoreType) { }
    void count_external_paired(size_t, size_t, ScoreType) { }
    void count_external_unpaired(size_t, size_t, ScoreType) { }
};
