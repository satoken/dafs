/*
 * LinearAlign implementation for DAFS
 */

#include "linearalign.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <unordered_map>

// Include the simplified BeamAlign implementation
#include "linearalign/BeamAlignSimple.h"

// Default HMM parameters based on common RNA alignment models
static const double DEFAULT_TRANS_PROBS[3][3] = {
    // M     X     Y
    {0.9,  0.05, 0.05},  // from M
    {0.15, 0.8,  0.05},  // from X  
    {0.15, 0.05, 0.8}    // from Y
};

static const double DEFAULT_EMIT_PROBS[3][5] = {
    // A     C     G     U     N
    {0.25, 0.25, 0.25, 0.25, 0.0},  // Match state
    {0.2,  0.2,  0.2,  0.2,  0.2},  // Insert X
    {0.2,  0.2,  0.2,  0.2,  0.2}   // Insert Y
};

LinearAlign::LinearAlign(float th, int beam_size)
    : Align::Model(th), 
      beam_size_(beam_size),
      use_prior_(false),
      trans_probs_(nullptr),
      emit_probs_(nullptr),
      custom_params_(false)
{
    // Create BeamAlign instance with specified beam size
    beam_align_ = std::make_unique<BeamAlign>(beam_size);
    
    // Initialize default parameters
    initializeDefaultParameters();
}

LinearAlign::~LinearAlign()
{
    // Clean up allocated HMM parameters if using defaults
    if (!custom_params_ && trans_probs_ != nullptr) {
        for (int i = 0; i < 3; ++i) {
            delete[] trans_probs_[i];
            delete[] emit_probs_[i];
        }
        delete[] trans_probs_;
        delete[] emit_probs_;
    }
}

void LinearAlign::initializeDefaultParameters()
{
    if (!custom_params_) {
        // Allocate and set default HMM parameters
        trans_probs_ = new double*[3];
        emit_probs_ = new double*[3];
        
        for (int i = 0; i < 3; ++i) {
            trans_probs_[i] = new double[3];
            emit_probs_[i] = new double[5];
            
            for (int j = 0; j < 3; ++j) {
                trans_probs_[i][j] = std::log(DEFAULT_TRANS_PROBS[i][j]);
            }
            for (int j = 0; j < 5; ++j) {
                emit_probs_[i][j] = std::log(DEFAULT_EMIT_PROBS[i][j]);
            }
        }
    }
}

void LinearAlign::setHMMParameters(double** trans_probs, double** emit_probs)
{
    // Clean up default parameters if they were allocated
    if (!custom_params_ && trans_probs_ != nullptr) {
        for (int i = 0; i < 3; ++i) {
            delete[] trans_probs_[i];
            delete[] emit_probs_[i];
        }
        delete[] trans_probs_;
        delete[] emit_probs_;
    }
    
    trans_probs_ = trans_probs;
    emit_probs_ = emit_probs;
    custom_params_ = true;
}

void LinearAlign::calculate(const std::string& seq1, const std::string& seq2, MP& mp)
{
    // Convert sequences to format expected by BeamAlign
    std::string seq1_copy = seq1;
    std::string seq2_copy = seq2;
    
    // Replace T with U for RNA sequences
    std::replace(seq1_copy.begin(), seq1_copy.end(), 'T', 'U');
    std::replace(seq2_copy.begin(), seq2_copy.end(), 'T', 'U');
    
    // Set up a simple match score function for BeamAlign
    // This can be improved with more sophisticated scoring
    beam_align_->setMatchScoreFunction([&seq1_copy, &seq2_copy](char a, char b) -> double {
        // Simple sequence similarity score
        if (a == b) {
            return 2.0;  // Match bonus
        } else if ((a == 'U' && b == 'T') || (a == 'T' && b == 'U')) {
            return 2.0;  // RNA T-U equivalence
        }
        return -1.0;  // Mismatch penalty
    });
    
    // Vectors to store alignment results
    std::vector<char> aln1, aln2;
    
    // Run BeamAlign ML alignment
    beam_align_->ml_alignment(seq1_copy, seq2_copy, aln1, aln2, 
                             trans_probs_, emit_probs_, use_prior_);
    
    // Run forward algorithm to get alignment score
    double forward_score = beam_align_->forward(seq1_copy, seq2_copy,
                                                trans_probs_, emit_probs_, use_prior_);
    
    // Run backward algorithm
    beam_align_->backward(trans_probs_, emit_probs_, use_prior_);
    
    // Calculate alignment probabilities with threshold
    std::unordered_map<int, aln_ret>* aln_results = nullptr;
    aln_results = beam_align_->cal_align_prob(forward_score, threshold(), aln_results);
    
    // Convert results to DAFS sparse matrix format
    if (aln_results != nullptr) {
        convertToSparseMatrix(aln_results, seq1, seq2, mp);
        delete aln_results;
    }
}

void LinearAlign::convertToSparseMatrix(const std::unordered_map<int, struct aln_ret>* aln_results,
                                       const std::string& seq1, const std::string& seq2,
                                       MP& mp)
{
    const uint L1 = seq1.size();
    const uint L2 = seq2.size();
    
    // Initialize sparse matrix
    mp.clear();
    mp.resize(L1);
    
    // Convert from BeamAlign's format to DAFS sparse matrix format
    // BeamAlign uses: key = i * max_len + j
    uint max_len = std::max(L1, L2) + 1;
    
    for (const auto& result : *aln_results) {
        int key = result.first;
        const struct aln_ret& aln_info = result.second;
        
        // Decode position from key
        uint i = key / max_len;
        uint j = key % max_len;
        
        // Check bounds and threshold
        if (i < L1 && j < L2) {
            // Convert from log probability to probability
            float prob = std::exp(aln_info.aln_prob);
            
            if (prob > threshold()) {
                mp[i].push_back(std::make_pair(j, prob));
            }
        }
    }
    
    // Sort each row by column index for consistency
    for (uint i = 0; i < L1; ++i) {
        std::sort(mp[i].begin(), mp[i].end(),
                 [](const std::pair<uint, float>& a, const std::pair<uint, float>& b) {
                     return a.first < b.first;
                 });
    }
}