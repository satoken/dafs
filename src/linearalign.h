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
#include "align.h"
#include "typedefs.h"

// Forward declarations
class BeamAlign;

class LinearAlign : public Align::Model
{
public:
    // Constructor with threshold and beam size
    LinearAlign(float th, int beam_size = 100);
    ~LinearAlign();
    
    // Main calculation method - implements Align::Model interface
    void calculate(const std::string& seq1, const std::string& seq2, MP& mp) override;
    
    // Set HMM parameters (transition and emission probabilities)
    void setHMMParameters(double** trans_probs, double** emit_probs);
    
    // Set whether to use prior information
    void setUsePrior(bool use_prior) { use_prior_ = use_prior; }
    
private:
    std::unique_ptr<BeamAlign> beam_align_;
    int beam_size_;
    bool use_prior_;
    
    // HMM parameters (if provided externally)
    double** trans_probs_;
    double** emit_probs_;
    bool custom_params_;
    
    // Helper method to convert BeamAlign output to DAFS sparse matrix format
    void convertToSparseMatrix(const std::unordered_map<int, struct aln_ret>* aln_results,
                               const std::string& seq1, const std::string& seq2,
                               MP& mp);
    
    // Initialize default HMM parameters if not provided
    void initializeDefaultParameters();
};

#endif // __INC_LINEARALIGN_H__