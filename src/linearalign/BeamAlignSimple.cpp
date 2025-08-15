/*
 * BeamAlignSimple.cpp
 * Simplified BeamAlign implementation for DAFS integration
 */

#include "BeamAlignSimple.h"
#include <algorithm>
#include <iostream>

BeamAlign::BeamAlign(int beam_size) 
    : beam(beam_size),
      match_score_func_(nullptr),
      match_score_(2.0),
      mismatch_score_(-1.0),
      gap_open_(-5.0),
      gap_extend_(-1.0)
{
}

BeamAlign::~BeamAlign() {
}

double BeamAlign::getMatchScore(char a, char b) {
    if (match_score_func_) {
        return match_score_func_(a, b);
    }
    
    // Default scoring
    if (a == b) {
        return match_score_;
    } else if ((a == 'U' && b == 'T') || (a == 'T' && b == 'U')) {
        return match_score_;  // RNA T-U equivalence
    } else {
        return mismatch_score_;
    }
}

double BeamAlign::getTransitionScore(int from_state, int to_state) {
    // Simple HMM transition scores
    // States: 0=match, 1=ins1, 2=ins2
    
    if (from_state == 0) {  // From match
        if (to_state == 0) return std::log(0.9);      // match->match
        if (to_state == 1) return std::log(0.05);     // match->ins1
        if (to_state == 2) return std::log(0.05);     // match->ins2
    } else if (from_state == 1) {  // From ins1
        if (to_state == 0) return std::log(0.1);      // ins1->match
        if (to_state == 1) return std::log(0.8);      // ins1->ins1
        if (to_state == 2) return std::log(0.1);      // ins1->ins2
    } else if (from_state == 2) {  // From ins2
        if (to_state == 0) return std::log(0.1);      // ins2->match
        if (to_state == 1) return std::log(0.1);      // ins2->ins1
        if (to_state == 2) return std::log(0.8);      // ins2->ins2
    }
    
    return std::log(1e-10);  // Very small probability
}

void BeamAlign::initializeMatrices(int len1, int len2) {
    forward_match_.assign(len1 + 1, std::vector<double>(len2 + 1, -1e9));
    forward_ins1_.assign(len1 + 1, std::vector<double>(len2 + 1, -1e9));
    forward_ins2_.assign(len1 + 1, std::vector<double>(len2 + 1, -1e9));
    
    backward_match_.assign(len1 + 1, std::vector<double>(len2 + 1, -1e9));
    backward_ins1_.assign(len1 + 1, std::vector<double>(len2 + 1, -1e9));
    backward_ins2_.assign(len1 + 1, std::vector<double>(len2 + 1, -1e9));
}

double BeamAlign::runForward() {
    int len1 = seq1_.length();
    int len2 = seq2_.length();
    
    // Initialize
    forward_match_[0][0] = 0.0;
    
    // Fill DP matrix
    for (int i = 0; i <= len1; ++i) {
        for (int j = 0; j <= len2; ++j) {
            if (i == 0 && j == 0) continue;
            
            // Match state
            if (i > 0 && j > 0) {
                double emit_score = getMatchScore(seq1_[i-1], seq2_[j-1]);
                forward_match_[i][j] = std::max({
                    forward_match_[i-1][j-1] + getTransitionScore(0, 0) + emit_score,
                    forward_ins1_[i-1][j-1] + getTransitionScore(1, 0) + emit_score,
                    forward_ins2_[i-1][j-1] + getTransitionScore(2, 0) + emit_score
                });
            }
            
            // Insert seq1 state
            if (i > 0) {
                double emit_score = gap_extend_;
                forward_ins1_[i][j] = std::max({
                    forward_match_[i-1][j] + getTransitionScore(0, 1) + emit_score,
                    forward_ins1_[i-1][j] + getTransitionScore(1, 1) + emit_score,
                    forward_ins2_[i-1][j] + getTransitionScore(2, 1) + emit_score
                });
            }
            
            // Insert seq2 state
            if (j > 0) {
                double emit_score = gap_extend_;
                forward_ins2_[i][j] = std::max({
                    forward_match_[i][j-1] + getTransitionScore(0, 2) + emit_score,
                    forward_ins1_[i][j-1] + getTransitionScore(1, 2) + emit_score,
                    forward_ins2_[i][j-1] + getTransitionScore(2, 2) + emit_score
                });
            }
        }
    }
    
    // Return final score
    return std::max({
        forward_match_[len1][len2],
        forward_ins1_[len1][len2],
        forward_ins2_[len1][len2]
    });
}

void BeamAlign::runBackward() {
    int len1 = seq1_.length();
    int len2 = seq2_.length();
    
    // Initialize final states
    backward_match_[len1][len2] = 0.0;
    backward_ins1_[len1][len2] = 0.0;
    backward_ins2_[len1][len2] = 0.0;
    
    // Fill backward matrix
    for (int i = len1; i >= 0; --i) {
        for (int j = len2; j >= 0; --j) {
            if (i == len1 && j == len2) continue;
            
            // Update from future states
            if (i < len1 && j < len2) {
                double emit_score = getMatchScore(seq1_[i], seq2_[j]);
                backward_match_[i][j] = std::max({
                    backward_match_[i+1][j+1] + getTransitionScore(0, 0) + emit_score,
                    backward_ins1_[i+1][j+1] + getTransitionScore(0, 1) + gap_extend_,
                    backward_ins2_[i+1][j+1] + getTransitionScore(0, 2) + gap_extend_
                });
            }
            
            if (i < len1) {
                backward_ins1_[i][j] = std::max({
                    backward_match_[i+1][j] + getTransitionScore(1, 0) + getMatchScore(seq1_[i], '-'),
                    backward_ins1_[i+1][j] + getTransitionScore(1, 1) + gap_extend_,
                    backward_ins2_[i+1][j] + getTransitionScore(1, 2) + gap_extend_
                });
            }
            
            if (j < len2) {
                backward_ins2_[i][j] = std::max({
                    backward_match_[i][j+1] + getTransitionScore(2, 0) + getMatchScore('-', seq2_[j]),
                    backward_ins1_[i][j+1] + getTransitionScore(2, 1) + gap_extend_,
                    backward_ins2_[i][j+1] + getTransitionScore(2, 2) + gap_extend_
                });
            }
        }
    }
}

void BeamAlign::computePosteriorProbs(std::unordered_map<int, aln_ret>* results, double threshold) {
    int len1 = seq1_.length();
    int len2 = seq2_.length();
    int max_len = std::max(len1, len2) + 1;
    
    double total_score = runForward();
    
    for (int i = 1; i <= len1; ++i) {
        for (int j = 1; j <= len2; ++j) {
            // Calculate posterior probability for match (i,j)
            double post_prob = forward_match_[i][j] + backward_match_[i][j] - total_score;
            
            if (post_prob > std::log(threshold)) {
                int key = i * max_len + j;
                aln_ret& ret = (*results)[key];
                ret.prob = post_prob;
                ret.aln_prob = post_prob;
            }
        }
    }
}

void BeamAlign::ml_alignment(std::string &seq1, std::string &seq2,
                            std::vector<char> &aln1, std::vector<char> &aln2,
                            double** &transprobs, double** &emitprobs, bool prior) {
    seq1_ = seq1;
    seq2_ = seq2;
    
    int len1 = seq1.length();
    int len2 = seq2.length();
    
    initializeMatrices(len1, len2);
    runForward();
    
    // Simple traceback for ML alignment
    aln1.clear();
    aln2.clear();
    
    int i = len1, j = len2;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0) {
            aln1.push_back(seq1[i-1]);
            aln2.push_back(seq2[j-1]);
            i--; j--;
        } else if (i > 0) {
            aln1.push_back(seq1[i-1]);
            aln2.push_back('-');
            i--;
        } else {
            aln1.push_back('-');
            aln2.push_back(seq2[j-1]);
            j--;
        }
    }
    
    std::reverse(aln1.begin(), aln1.end());
    std::reverse(aln2.begin(), aln2.end());
}

double BeamAlign::forward(std::string seq1, std::string seq2,
                         double** &trans_probs, double** &emit_probs, bool prior) {
    seq1_ = seq1;
    seq2_ = seq2;
    
    initializeMatrices(seq1.length(), seq2.length());
    return runForward();
}

double BeamAlign::backward(double** &transprobs, double** &emitprobs, bool prior) {
    runBackward();
    return 0.0;
}

std::unordered_map<int, aln_ret>* BeamAlign::cal_align_prob(double forward_score,
                                                           double threshold,
                                                           std::unordered_map<int, aln_ret>* &aln_results) {
    auto* results = new std::unordered_map<int, aln_ret>;
    computePosteriorProbs(results, threshold);
    return results;
}