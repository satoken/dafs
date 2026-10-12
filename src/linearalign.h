/*
 * LinearAlign wrapper for DAFS
 * 
 * This file integrates the BeamAlign algorithm from LinearTurboFold
 * into the DAFS alignment framework.
 */

#ifndef __INC_LINEARALIGN_H__
#define __INC_LINEARALIGN_H__

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>
#include <array>
#include "align.h"
#include "typedefs.h"

// Forward declarations
class BeamAlign;

class LinearAlign : public Align::Model
{
public:
    // CONTRAlign keeps its complete five-state, double-affine RNA model.
    // The other parameter sets use BeamAlign's three states.
    enum class ScoreModel {
        LinearTurboFold,
        CONTRAlign,
        ProbConsRNA
    };

    // Constructor with threshold and beam size
    LinearAlign(float th, int beam_size = 100,
                ScoreModel score_model = ScoreModel::LinearTurboFold);
    ~LinearAlign();

    static ScoreModel parseScoreModel(const std::string& name);
    static const char* scoreModelName(ScoreModel model);
    ScoreModel scoreModel() const { return score_model_; }
    unsigned scoreStateCount() const { return score_model_ == ScoreModel::CONTRAlign ? 5 : 3; }
    
    // Main calculation method - implements Align::Model interface
    void calculate(const std::string& seq1, const std::string& seq2, MP& mp) override;
    
    // Set HMM parameters (transition and emission probabilities)
    void setHMMParameters(double** trans_probs, double** emit_probs);
    
    // Set whether to use prior information
    void setUsePrior(bool use_prior) { use_prior_ = use_prior; }
    
private:
    using TransitionScores = std::array<std::array<double, 3>, 3>;
    using EmissionScores = std::array<std::array<double, 3>, 27>;

    std::unique_ptr<BeamAlign> beam_align_;
    int beam_size_;
    bool use_prior_;
    ScoreModel score_model_;
    TransitionScores transition_scores_;
    EmissionScores emission_scores_;
    bool parameters_initialized_;
    
    // Helper method to convert BeamAlign output to DAFS sparse matrix format
    void convertToSparseMatrix(const std::unordered_map<int, struct aln_ret>* aln_results,
                               const std::string& seq1, const std::string& seq2,
                               MP& mp);
    
    // Initialize the selected scoring model without reducing its state space.
    void initializeParameters(ScoreModel model);
};

#endif // __INC_LINEARALIGN_H__
