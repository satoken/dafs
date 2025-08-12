#pragma once

#include <string>
#include <vector>
#include "nested_array.hpp"

class CONTRAfoldNearestNeighbor
{
    public:
        using ScoreType = float;

    private:
        using SeqType = std::vector<short>;

    public:
        CONTRAfoldNearestNeighbor(const std::string& seq);
        ~CONTRAfoldNearestNeighbor() {};

        ScoreType score_hairpin(size_t i, size_t j) const;
        ScoreType score_single_loop(size_t i, size_t j, size_t k, size_t l) const;
        ScoreType score_helix(size_t i, size_t j, size_t m) const;
        ScoreType score_multi_loop(size_t i, size_t j) const;
        ScoreType score_multi_paired(size_t i, size_t j) const;
        ScoreType score_multi_unpaired(size_t i, size_t j) const;
        ScoreType score_external_zero() const { return 0.0; }
        ScoreType score_external_paired(size_t i, size_t j) const;
        ScoreType score_external_unpaired(size_t i, size_t j) const;

        void count_hairpin(size_t i, size_t j, ScoreType v);
        void count_single_loop(size_t i, size_t j, size_t k, size_t l, ScoreType v);
        void count_helix(size_t i, size_t j, size_t m, ScoreType v);
        void count_multi_loop(size_t i, size_t j, ScoreType v);
        void count_multi_paired(size_t i, size_t j, ScoreType v);
        void count_multi_unpaired(size_t i, size_t j, ScoreType v);
        void count_external_zero(ScoreType v) { }
        void count_external_paired(size_t i, size_t j, ScoreType v);
        void count_external_unpaired(size_t i, size_t j, ScoreType v);

    private:
        static auto convert_sequence(const std::string& seq) -> SeqType;

#if 0
        ScoreType score_base_pair(short i, short j) const;
        ScoreType score_helix_stacking(short i, short j, short k, short l) const;
        ScoreType score_internal_1x1(short i, short j) const;
#endif
        void cache_count_base_pair(short i, short j, ScoreType v);
        void cache_count_helix_stacking(short i, short j, short k, short l, ScoreType v);
        void cache_count_internal_1x1(short i, short j, ScoreType v);

    private:
        SeqType seq2_;

        array_2d<float, 5, 5> score_base_pair_;
        array_4d<float, 5, 5, 5, 5> score_terminal_mismatch_;
        array_1d<float, 31> score_hairpin_length_;
        array_2d<float, 5, 5> score_internal_explicit_;
        array_1d<float, 31> score_bulge_length_;
        array_1d<float, 31> score_internal_length_;
        array_1d<float, 16> score_internal_symmetry_;
        array_1d<float, 29> score_internal_asymmetry_;
        array_1d<float, 5> score_bulge_0x1_;
        array_2d<float, 5, 5> score_internal_1x1_;
        array_4d<float, 5, 5, 5, 5> score_helix_stacking_;
        array_2d<float, 5, 5> score_helix_closing_;
        array_1d<float, 1> score_multi_base_;
        array_1d<float, 1> score_multi_unpaired_;
        array_1d<float, 1> score_multi_paired_;
        array_3d<float, 5, 5, 5> score_dangle_left_;
        array_3d<float, 5, 5, 5> score_dangle_right_;
        array_1d<float, 1> score_external_unpaired_;
        array_1d<float, 1> score_external_paired_;

        array_2d<float, 5, 5> count_base_pair_;
        array_4d<float, 5, 5, 5, 5> count_terminal_mismatch_;
        array_1d<float, 31> count_hairpin_length_;
        array_2d<float, 5, 5> count_internal_explicit_;
        array_1d<float, 31> count_bulge_length_;
        array_1d<float, 31> count_internal_length_;
        array_1d<float, 16> count_internal_symmetry_;
        array_1d<float, 29> count_internal_asymmetry_;
        array_1d<float, 5> count_bulge_0x1_;
        array_2d<float, 5, 5> count_internal_1x1_;
        array_4d<float, 5, 5, 5, 5> count_helix_stacking_;
        array_2d<float, 5, 5> count_helix_closing_;
        array_1d<float, 1> count_multi_base_;
        array_1d<float, 1> count_multi_unpaired_;
        array_1d<float, 1> count_multi_paired_;
        array_3d<float, 5, 5, 5> count_dangle_left_;
        array_3d<float, 5, 5, 5> count_dangle_right_;
        array_1d<float, 1> count_external_unpaired_;
        array_1d<float, 1> count_external_paired_;

        std::vector<float> cache_score_hairpin_length_;
        std::vector<float> cache_score_bulge_length_;
        std::vector<float> cache_score_internal_length_;
        std::vector<float> cache_score_internal_symmetry_;
        std::vector<float> cache_score_internal_asymmetry_;
        std::vector<std::vector<float>> cache_score_base_pair_;
        std::vector<std::vector<std::vector<std::vector<float>>>> cache_score_helix_stacking_;
        std::vector<std::vector<float>> cache_score_internal_1x1_;

    public:
        const u_int32_t MAX_HAIRPIN_LENGTH;
        const u_int32_t MAX_BULGE_LENGTH;
        const u_int32_t MAX_INTERNAL_LENGTH;
        const u_int32_t MAX_SINGLE_LENGTH;
        const u_int32_t MAX_INTERNAL_SYMMETRIC_LENGTH;
        const u_int32_t MAX_INTERNAL_ASYMMETRY;
        const u_int32_t MAX_INTERNAL_EXPLICIT_LENGTH;
};
