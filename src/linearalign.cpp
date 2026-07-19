/*
 * LinearAlign implementation for DAFS
 * Uses BeamAlign with ML parameters from LinearTurboFold for posterior probability calculation
 */

#include "linearalign.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <unordered_map>
#include <iomanip>
#include <stdexcept>

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

LinearAlign::LinearAlign(float th, int beam_size)
    : Align::Model(th), 
      beam_size_(beam_size),
      use_prior_(false),  // Use false for initial iteration (no structure info)
      trans_probs_(nullptr),
      emit_probs_(nullptr),
      custom_params_(false)
{
    beam_align_ = std::make_unique<BeamAlign>(beam_size_);
    initializeMLParameters();
}

LinearAlign::~LinearAlign()
{
    cleanupParameters();
}

void LinearAlign::initializeMLParameters()
{
    // Clean up any existing parameters
    cleanupParameters();
    
    // Allocate and initialize ML transition probabilities
    trans_probs_ = new double*[3];
    for (int i = 0; i < 3; ++i) {
        trans_probs_[i] = new double[3];
        for (int j = 0; j < 3; ++j) {
            trans_probs_[i][j] = std::log(ML_TRANS_PROBS[i][j]);
        }
    }
    
    // Allocate and initialize ML emission probabilities
    emit_probs_ = new double*[27];
    for (int sym = 0; sym < 27; ++sym) {
        emit_probs_[sym] = new double[3];
        for (int state = 0; state < 3; ++state) {
            emit_probs_[sym][state] = std::log(ML_EMIT_PROBS[sym][state]);
        }
    }
    
    custom_params_ = true;
}

void LinearAlign::cleanupParameters()
{
    if (custom_params_) {
        if (trans_probs_) {
            for (int i = 0; i < 3; ++i) {
                delete[] trans_probs_[i];
            }
            delete[] trans_probs_;
            trans_probs_ = nullptr;
        }
        
        if (emit_probs_) {
            for (int i = 0; i < 27; ++i) {
                delete[] emit_probs_[i];
            }
            delete[] emit_probs_;
            emit_probs_ = nullptr;
        }
        
        custom_params_ = false;
    }
}

void LinearAlign::setHMMParameters(double** trans_probs, double** emit_probs)
{
    // Clean up default parameters if they were allocated
    if (custom_params_) {
        cleanupParameters();
    }
    
    trans_probs_ = trans_probs;
    emit_probs_ = emit_probs;
    custom_params_ = false;  // External parameters, don't clean up in destructor
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
        if (!trans_probs_ || !emit_probs_)
            throw std::logic_error("HMM parameters are not initialized");

        // Step 1: Run forward algorithm
        double forward_score = beam_align_->forward(seq1_copy, seq2_copy, trans_probs_, emit_probs_, use_prior_);
        if (!std::isfinite(forward_score))
            throw std::runtime_error("forward algorithm returned a non-finite score");
        
        // Step 2: Run backward algorithm  
        double backward_score = beam_align_->backward(trans_probs_, emit_probs_, use_prior_);
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
            "LinearAlign failed for sequence lengths " +
            std::to_string(seq1.size()) + " and " +
            std::to_string(seq2.size()) + ": " + e.what());
    } catch (...) {
        throw std::runtime_error(
            "LinearAlign failed for sequence lengths " +
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

// Keep the old initializeDefaultParameters method for compatibility
void LinearAlign::initializeDefaultParameters()
{
    // This method is now replaced by initializeMLParameters
    // but kept for compatibility
    initializeMLParameters();
}
