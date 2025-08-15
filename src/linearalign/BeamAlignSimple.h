/*
 * BeamAlignSimple.h
 * Simplified BeamAlign implementation for DAFS integration
 */

#ifndef BEAM_ALIGN_SIMPLE_H
#define BEAM_ALIGN_SIMPLE_H

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <cmath>

// Simplified structures for alignment
struct aln_ret {
    double prob;
    double aln_prob;
    
    aln_ret() : prob(0.0), aln_prob(0.0) {}
};

struct AlignState {
    double score;
    int manner;  // 0=match, 1=ins1, 2=ins2
    int prev_i, prev_j;
    
    AlignState() : score(-1e9), manner(-1), prev_i(-1), prev_j(-1) {}
    AlignState(double s, int m, int pi, int pj) : score(s), manner(m), prev_i(pi), prev_j(pj) {}
};

// Score function interface
typedef std::function<double(char, char)> MatchScoreFunction;

class BeamAlign {
public:
    int beam;
    
    BeamAlign(int beam_size = 100);
    ~BeamAlign();
    
    // Set custom match scoring function
    void setMatchScoreFunction(MatchScoreFunction func) { match_score_func_ = func; }
    
    // Main alignment methods
    void ml_alignment(std::string &seq1, std::string &seq2, 
                     std::vector<char> &aln1, std::vector<char> &aln2,
                     double** &transprobs, double** &emitprobs, bool prior);
    
    double forward(std::string seq1, std::string seq2,
                  double** &trans_probs, double** &emit_probs, bool prior);
    
    double backward(double** &transprobs, double** &emitprobs, bool prior);
    
    std::unordered_map<int, aln_ret>* cal_align_prob(double forward_score, 
                                                     double threshold,
                                                     std::unordered_map<int, aln_ret>* &aln_results);

private:
    MatchScoreFunction match_score_func_;
    
    // Cached sequences for alignment
    std::string seq1_, seq2_;
    
    // DP matrices for forward/backward
    std::vector<std::vector<double>> forward_match_, forward_ins1_, forward_ins2_;
    std::vector<std::vector<double>> backward_match_, backward_ins1_, backward_ins2_;
    
    // Default scoring parameters
    double match_score_;
    double mismatch_score_;
    double gap_open_;
    double gap_extend_;
    
    // Helper methods
    double getMatchScore(char a, char b);
    double getTransitionScore(int from_state, int to_state);
    void initializeMatrices(int len1, int len2);
    
    // DP algorithms
    double runForward();
    void runBackward();
    void computePosteriorProbs(std::unordered_map<int, aln_ret>* results, double threshold);
};

#endif // BEAM_ALIGN_SIMPLE_H