/*
 * LinearAlign implementation for DAFS
 * Uses BeamAlign with selectable LinearTurboFold, CONTRAlign, or ProbConsRNA
 * scoring parameters for posterior probability calculation.
 */

#include "linearalign.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <stdexcept>
#include <cctype>
#include <limits>

// Include the original BeamAlign implementation
#include "linearalign/BeamAlign.h"

// ML HMM parameters from LinearTurboFold (trained on RNA families)
// States: 0=INS1, 1=INS2, 2=ALIGN
static const double ML_TRANS_PROBS[3][3] = {
    {0.666439, 0.041319, 0.292242}, // INS1 -> [INS1, INS2, ALIGN]
    {0.041319, 0.666439, 0.292242}, // INS2 -> [INS1, INS2, ALIGN]
    {0.022666, 0.022666, 0.954668}  // ALIGN -> [INS1, INS2, ALIGN]
};

// ML emission probabilities from LinearTurboFold
// 27 symbols: 25 nucleotide pairs (5x5) + start(25) + end(26)
// Order: AA, AC, AG, AU, A., CA, CC, CG, CU, C., GA, GC, GG, GU, G., UA, UC, UG, UU, U., .A, .C, .G, .U, .., START, END
static const double ML_EMIT_PROBS[27][3] = {
    {0.000000, 0.000000, 0.134009}, // AA
    {0.000000, 0.000000, 0.027164}, // AC
    {0.000000, 0.000000, 0.049659}, // AG
    {0.000000, 0.000000, 0.028825}, // AU
    {0.211509, 0.000000, 0.000000}, // A.
    {0.000000, 0.000000, 0.027164}, // CA
    {0.000000, 0.000000, 0.140242}, // CC
    {0.000000, 0.000000, 0.037862}, // CG
    {0.000000, 0.000000, 0.047735}, // CU
    {0.257349, 0.000000, 0.000000}, // C.
    {0.000000, 0.000000, 0.049659}, // GA
    {0.000000, 0.000000, 0.037862}, // GC
    {0.000000, 0.000000, 0.178863}, // GG
    {0.000000, 0.000000, 0.032351}, // GU
    {0.271398, 0.000000, 0.000000}, // G.
    {0.000000, 0.000000, 0.028825}, // UA
    {0.000000, 0.000000, 0.047735}, // UC
    {0.000000, 0.000000, 0.032351}, // UG
    {0.000000, 0.000000, 0.099694}, // UU
    {0.259744, 0.000000, 0.000000}, // U.
    {0.000000, 0.211509, 0.000000}, // .A
    {0.000000, 0.257349, 0.000000}, // .C
    {0.000000, 0.271398, 0.000000}, // .G
    {0.000000, 0.259744, 0.000000}, // .U
    {0.000000, 0.000000, 0.000000}, // ..
    {0.000000, 0.000000, 1.000000}, // START
    {0.000000, 0.000000, 1.000000}  // END
};

namespace
{
constexpr double CONTRALIGN_MATCH[4][4] = {
    { 0.5256508867, -0.4090640200, -0.2502759109, -0.3252306723 },
    {-0.4090640200,  0.6665219366, -0.3289391181, -0.1326088918 },
    {-0.2502759109, -0.3289391181,  0.6684676551, -0.3565888168 },
    {-0.3252306723, -0.1326088918, -0.3565888168,  0.4590520450 }
};

constexpr double CONTRALIGN_INSERT[4] = {
    -0.0025219272, -0.0831389156, -0.0744397065, -0.0129005460
};

constexpr double PROBCONS_SINGLE[4] = {
    0.2270790040, 0.2422080040, 0.2839320004, 0.2464679927
};

constexpr double PROBCONS_PAIR[4][4] = {
    {0.1487240046, 0.0184142999, 0.0361397006, 0.0238473993},
    {0.0184142999, 0.1583919972, 0.0275536999, 0.0389291011},
    {0.0361397006, 0.0275536999, 0.1979320049, 0.0244289003},
    {0.0238473993, 0.0389291011, 0.0244289003, 0.1557479948}
};

double logScore(double value)
{
    return value > 0.0 ? std::log(value)
                       : -std::numeric_limits<double>::infinity();
}

std::string normalizedModelName(const std::string& name)
{
    std::string normalized;
    normalized.reserve(name.size());
    for (const unsigned char ch : name)
        if (std::isalnum(ch))
            normalized.push_back(static_cast<char>(std::tolower(ch)));
    return normalized;
}
} // namespace

LinearAlign::LinearAlign(float th, int beam_size, ScoreModel score_model)
    : Align::Model(th), 
      beam_size_(beam_size),
      use_prior_(false),  // Use false for initial iteration (no structure info)
      score_model_(score_model),
      parameters_initialized_(false)
{
    beam_align_ = std::make_unique<BeamAlign>(beam_size_);
    initializeParameters(score_model_);
}

LinearAlign::~LinearAlign() = default;

LinearAlign::ScoreModel LinearAlign::parseScoreModel(const std::string& name)
{
    const std::string normalized = normalizedModelName(name);
    if (normalized == "linearturbofold" || normalized == "turbofold")
        return ScoreModel::LinearTurboFold;
    if (normalized == "contralign")
        return ScoreModel::CONTRAlign;
    if (normalized == "probconsrna" || normalized == "probcons")
        return ScoreModel::ProbConsRNA;
    throw std::invalid_argument(
        "unknown LinearAlign score model '" + name +
        "' (expected LinearTurboFold, CONTRAlign, or ProbConsRNA)");
}

const char* LinearAlign::scoreModelName(ScoreModel model)
{
    switch (model) {
    case ScoreModel::LinearTurboFold: return "LinearTurboFold";
    case ScoreModel::CONTRAlign: return "CONTRAlign";
    case ScoreModel::ProbConsRNA: return "ProbConsRNA";
    }
    throw std::logic_error("invalid LinearAlign score model");
}

void LinearAlign::initializeParameters(ScoreModel model)
{
    // BeamAlign's xlog_sum implementation does not treat (-inf, -inf) as a
    // special case.  Use the same negligible floor as the bundled
    // ProbConsRNA wrapper for forbidden transitions and emissions.
    const double impossible = logScore(1e-10);
    for (auto& row : transition_scores_)
        row.fill(impossible);
    for (auto& row : emission_scores_)
        row.fill(impossible);

    if (model == ScoreModel::LinearTurboFold) {
        for (size_t from = 0; from < 3; ++from)
            for (size_t to = 0; to < 3; ++to)
                transition_scores_[from][to] = logScore(ML_TRANS_PROBS[from][to]);
        for (size_t symbol = 0; symbol < 27; ++symbol)
            for (size_t state = 0; state < 3; ++state)
                emission_scores_[symbol][state] = logScore(ML_EMIT_PROBS[symbol][state]);
    } else if (model == ScoreModel::ProbConsRNA) {
        constexpr double gap_open = 0.0190259293;
        constexpr double gap_extend = 0.3269913495;
        transition_scores_ = {{
            {{logScore(gap_extend), impossible, logScore(1.0 - gap_extend)}},
            {{impossible, logScore(gap_extend), logScore(1.0 - gap_extend)}},
            {{logScore(gap_open), logScore(gap_open), logScore(1.0 - 2.0 * gap_open)}}
        }};
        for (size_t a = 0; a < 4; ++a) {
            emission_scores_[a * 5 + 4][0] = logScore(PROBCONS_SINGLE[a]);
            emission_scores_[4 * 5 + a][1] = logScore(PROBCONS_SINGLE[a]);
            for (size_t b = 0; b < 4; ++b)
                emission_scores_[a * 5 + b][2] = logScore(PROBCONS_PAIR[a][b]);
        }
        emission_scores_[25][2] = 0.0;
        emission_scores_[26][2] = 0.0;
    } else {
        // CONTRAlign is a log-linear model.  BeamAlign consumes additive
        // log scores, so its single-affine subset maps directly onto the
        // three BeamAlign states.  CONTRAlign's second affine gap pair is
        // intentionally omitted because BeamAlign has exactly three states.
        constexpr double match_state = 0.3959924457;
        constexpr double insert_state = -0.4431756229;
        constexpr double match_to_match = 2.5057567100;
        constexpr double match_to_insert = -1.2423961130;
        constexpr double insert_extend = 1.8676346730;
        constexpr double insert_change = -6.9696754440;
        transition_scores_ = {{
            {{insert_extend, insert_change, match_to_insert}},
            {{insert_change, insert_extend, match_to_insert}},
            {{match_to_insert, match_to_insert, match_to_match}}
        }};
        for (size_t a = 0; a < 4; ++a) {
            emission_scores_[a * 5 + 4][0] = CONTRALIGN_INSERT[a] + insert_state;
            emission_scores_[4 * 5 + a][1] = CONTRALIGN_INSERT[a] + insert_state;
            for (size_t b = 0; b < 4; ++b)
                emission_scores_[a * 5 + b][2] = CONTRALIGN_MATCH[a][b] + match_state;
        }
        emission_scores_[25][2] = 0.0;
        emission_scores_[26][2] = 0.0;
    }

    score_model_ = model;
    parameters_initialized_ = true;
}

void LinearAlign::setHMMParameters(double** trans_probs, double** emit_probs)
{
    if (!trans_probs || !emit_probs) {
        parameters_initialized_ = false;
        return;
    }
    for (size_t from = 0; from < 3; ++from)
        for (size_t to = 0; to < 3; ++to)
            transition_scores_[from][to] = trans_probs[from][to];
    for (size_t symbol = 0; symbol < 27; ++symbol)
        for (size_t state = 0; state < 3; ++state)
            emission_scores_[symbol][state] = emit_probs[symbol][state];
    parameters_initialized_ = true;
}

void LinearAlign::calculate(const std::string& seq1, const std::string& seq2, MP& mp)
{
    // Initialize empty matrix first
    mp.clear();
    mp.resize(seq1.length());
    
    // Check for empty sequences
    if (seq1.empty() || seq2.empty()) {
        return;
    }
    
    // Convert sequences to format expected by BeamAlign
    std::string seq1_copy = seq1;
    std::string seq2_copy = seq2;
    
    // Replace T with U for RNA sequences
    std::replace(seq1_copy.begin(), seq1_copy.end(), 'T', 'U');
    std::replace(seq2_copy.begin(), seq2_copy.end(), 'T', 'U');
    
    try {
        // Invalid model state is a hard failure.  Continuing with an empty
        // posterior matrix makes downstream output look successful.
        if (!parameters_initialized_)
            throw std::logic_error("HMM parameters are not initialized");

        std::array<double*, 3> transition_rows;
        std::array<double*, 27> emission_rows;
        for (size_t i = 0; i < transition_rows.size(); ++i)
            transition_rows[i] = transition_scores_[i].data();
        for (size_t i = 0; i < emission_rows.size(); ++i)
            emission_rows[i] = emission_scores_[i].data();
        double** trans_probs = transition_rows.data();
        double** emit_probs = emission_rows.data();

        // Step 1: Run forward algorithm
        double forward_score = beam_align_->forward(seq1_copy, seq2_copy, trans_probs, emit_probs, use_prior_);
        if (!std::isfinite(forward_score))
            throw std::runtime_error("forward algorithm returned a non-finite score");
        
        // Step 2: Run backward algorithm  
        double backward_score = beam_align_->backward(trans_probs, emit_probs, use_prior_);
        if (!std::isfinite(backward_score))
            throw std::runtime_error("backward algorithm returned a non-finite score");
        
        // Step 3: Calculate posterior alignment probabilities
        std::unordered_map<int, aln_ret>* aln_results = nullptr;  // Always start with nullptr
        double log_threshold = std::log(this->threshold());  // Convert threshold to log space
        
        aln_results = beam_align_->cal_align_prob(forward_score, log_threshold, aln_results);
        
        // Step 4: Convert results to DAFS sparse matrix format
        if (aln_results == nullptr)
            throw std::runtime_error("posterior calculation returned no result matrix");
        std::unique_ptr<std::unordered_map<int, aln_ret>[]> result_owner(aln_results);
        convertToSparseMatrix(result_owner.get(), seq1_copy, seq2_copy, mp);
        
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "LinearAlign failed (score model " +
            std::string(scoreModelName(score_model_)) + ") for sequence lengths " +
            std::to_string(seq1.size()) + " and " +
            std::to_string(seq2.size()) + ": " + e.what());
    } catch (...) {
        throw std::runtime_error(
            "LinearAlign failed (score model " +
            std::string(scoreModelName(score_model_)) + ") for sequence lengths " +
            std::to_string(seq1.size()) + " and " +
            std::to_string(seq2.size()) + ": unknown exception");
    }
}

void LinearAlign::convertToSparseMatrix(const std::unordered_map<int, aln_ret>* aln_results,
                                       const std::string& seq1, const std::string& seq2,
                                       MP& mp)
{
    const uint L1 = seq1.size();
    const uint L2 = seq2.size();
    
    // Initialize sparse matrix
    mp.clear();
    mp.resize(L1);
    
    // Convert from BeamAlign's format to DAFS sparse matrix format
    // BeamAlign uses seq_len = original_len + 1 for start/end symbols
    // aln_results array size = seq1_len = L1 + 1
    // Valid array indices: 0 to L1 (inclusive)
    // However, BeamAlign typically uses indices 1 to L1 for actual sequence positions
    
    uint seq1_len = L1 + 1;  // BeamAlign's internal sequence length
    
    // Only iterate through valid sequence positions (1-based in BeamAlign)
    for (uint i = 1; i <= L1; ++i) {
        // Ensure we don't access beyond allocated array size
        if (i >= seq1_len) break;
        
        const auto& pos_results = aln_results[i];
        for (const auto& result : pos_results) {
            uint j = result.first;  // j is position in seq2 (1-based in BeamAlign)
            const aln_ret& aln_info = result.second;
            
            // Validate j is within valid sequence range (1-based)
            if (j >= 1 && j <= L2) {
                float prob = aln_info.prob;
                
                // Validate probability value (must be finite and reasonable)
                if (std::isfinite(prob) && !std::isnan(prob) && prob > 0.0f && prob <= 1.0f && prob > this->threshold()) {
                    // Convert to 0-based indexing for DAFS
                    uint i_zero = i - 1;  // Convert from 1-based to 0-based
                    uint j_zero = j - 1;  // Convert from 1-based to 0-based
                    
                    // Final bounds check before adding to sparse matrix
                    if (i_zero < L1 && j_zero < L2 && i_zero < mp.size()) {
                        mp[i_zero].push_back(std::make_pair(j_zero, prob));
                    }
                }
            }
        }
    }
    
    // Sort each row by column index for DAFS compatibility
    for (uint i = 0; i < mp.size(); ++i) {
        std::sort(mp[i].begin(), mp[i].end(),
                 [](const std::pair<uint, float>& a, const std::pair<uint, float>& b) {
                     return a.first < b.first;
                 });
    }
}
