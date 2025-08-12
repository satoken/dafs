#pragma once

#include <string>
#include <vector>
#include "nested_array.hpp"

class TurnerNearestNeighbor
{
    public:
        using ScoreType = float;

    private:
        using SeqType = std::vector<short>;

    public:
        TurnerNearestNeighbor(const std::string& seq);
        ~TurnerNearestNeighbor() {};

        ScoreType score_hairpin(size_t i, size_t j) const;
        ScoreType score_single_loop(size_t i, size_t j, size_t k, size_t l) const;
        ScoreType score_helix(size_t i, size_t j, size_t m) const;
        ScoreType score_multi_loop(size_t i, size_t j) const;
        ScoreType score_multi_paired(size_t i, size_t j) const;
        ScoreType score_multi_unpaired(size_t i, size_t j) const;
        ScoreType score_external_zero() const { return 0.0; }
        ScoreType score_external_paired(size_t i, size_t j) const;
        ScoreType score_external_unpaired(size_t i, size_t j) const { return 0.0; }

        void count_hairpin(size_t i, size_t j, ScoreType v);
        void count_single_loop(size_t i, size_t j, size_t k, size_t l, ScoreType v);
        void count_helix(size_t i, size_t j, size_t m, ScoreType v);
        void count_multi_loop(size_t i, size_t j, ScoreType v);
        void count_multi_paired(size_t i, size_t j, ScoreType v);
        void count_multi_unpaired(size_t i, size_t j, ScoreType v);
        void count_external_zero(ScoreType v) { }
        void count_external_paired(size_t i, size_t j, ScoreType v);
        void count_external_unpaired(size_t i, size_t j, ScoreType v) { }

    private:
        static auto convert_sequence(const std::string& seq) -> SeqType;

    private:
        SeqType seq2_;

        bool use_score_hairpin_at_least_;
        bool use_score_bulge_at_least_;
        bool use_score_internal_at_least_;

        array_2d<float, 8, 8> score_stack_;
        array_1d<float, 31> score_hairpin_;
        array_1d<float, 31> score_bulge_;
        array_1d<float, 31> score_internal_;
        array_3d<float, 8, 5, 5> score_mismatch_external_;
        array_3d<float, 8, 5, 5> score_mismatch_hairpin_;
        array_3d<float, 8, 5, 5> score_mismatch_internal_;
        array_3d<float, 8, 5, 5> score_mismatch_internal_1n_;
        array_3d<float, 8, 5, 5> score_mismatch_internal_23_;
        array_3d<float, 8, 5, 5> score_mismatch_multi_;
        array_4d<float, 8, 8, 5, 5> score_int11_;
        array_5d<float, 8, 8, 5, 5, 5> score_int21_;
        array_6d<float, 7, 7, 5, 5, 5, 5> score_int22_;
        array_2d<float, 8, 5> score_dangle5_;
        array_2d<float, 8, 5> score_dangle3_;
        array_1d<float, 1> score_ml_base_;
        array_1d<float, 1> score_ml_closing_;
        array_1d<float, 1> score_ml_intern_;
        array_1d<float, 1> score_ninio_;
        array_1d<float, 1> score_max_ninio_;
        array_1d<float, 1> score_duplex_init_;
        array_1d<float, 1> score_terminalAU_;
        array_1d<float, 1> score_lxc_;

        bool use_count_hairpin_at_least_;
        bool use_count_bulge_at_least_;
        bool use_count_internal_at_least_;
        array_2d<float, 8, 8> count_stack_;
        array_1d<float, 31> count_hairpin_;
        array_1d<float, 31> count_bulge_;
        array_1d<float, 31> count_internal_;
        array_3d<float, 8, 5, 5> count_mismatch_external_;
        array_3d<float, 8, 5, 5> count_mismatch_hairpin_;
        array_3d<float, 8, 5, 5> count_mismatch_internal_;
        array_3d<float, 8, 5, 5> count_mismatch_internal_1n_;
        array_3d<float, 8, 5, 5> count_mismatch_internal_23_;
        array_3d<float, 8, 5, 5> count_mismatch_multi_;
        array_4d<float, 8, 8, 5, 5> count_int11_;
        array_5d<float, 8, 8, 5, 5, 5> count_int21_;
        array_6d<float, 7, 7, 5, 5, 5, 5> count_int22_;
        array_2d<float, 8, 5> count_dangle5_;
        array_2d<float, 8, 5> count_dangle3_;
        array_1d<float, 1> count_ml_base_;
        array_1d<float, 1> count_ml_closing_;
        array_1d<float, 1> count_ml_intern_;
        array_1d<float, 1> count_ninio_;
        array_1d<float, 1> count_max_ninio_;
        array_1d<float, 1> count_duplex_init_;
        array_1d<float, 1> count_terminalAU_;
        array_1d<float, 1> count_lxc_;

        std::vector<float> cache_score_hairpin_;
        std::vector<float> cache_score_bulge_;
        std::vector<float> cache_score_internal_;

    private:
        static int complement_pair[5][5];
};
