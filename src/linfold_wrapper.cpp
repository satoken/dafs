/*
 * LinFold wrapper implementation for DAFS
 */

#include "linfold_wrapper.h"
#include <algorithm>
#include <cassert>

LinFoldWrapper::LinFoldWrapper(float th, ModelType model_type, uint32_t beam_size)
  : Fold::Model(th), beam_size_(beam_size), model_type_(model_type), 
    linfold_turner_(nullptr), linfold_contra_(nullptr)
{
}

void 
LinFoldWrapper::calculate(const std::string& seq, BP& bp)
{
  try {
    std::vector<std::vector<std::pair<u_int32_t, float>>> bpp;
    
    if (model_type_ == ModelType::LPV) {
      // Initialize LinFold with TurnerNearestNeighbor parameters for this sequence
      auto turner_params = std::make_unique<TurnerNearestNeighbor>(seq);
      linfold_turner_ = std::make_unique<LinFold<TurnerNearestNeighbor>>(std::move(turner_params));
      
      // Set up options with beam size
      LinFold<TurnerNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      
      // First compute inside algorithm (forward pass)
      linfold_turner_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_turner_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities
      bpp = linfold_turner_->compute_basepairing_probabilities(seq, opt);
    } else { // ModelType::LPC
      // Initialize LinFold with CONTRAfoldNearestNeighbor parameters for this sequence
      auto contra_params = std::make_unique<CONTRAfoldNearestNeighbor>(seq);
      linfold_contra_ = std::make_unique<LinFold<CONTRAfoldNearestNeighbor>>(std::move(contra_params));
      
      // Set up options with beam size
      LinFold<CONTRAfoldNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      
      // First compute inside algorithm (forward pass)
      linfold_contra_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_contra_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities
      bpp = linfold_contra_->compute_basepairing_probabilities(seq, opt);
    }
    
    // Convert LinFold output format to BP format used by DAFS
    uint L = seq.size();
    bp.clear();
    bp.resize(L);
    
    for (uint i = 0; i < L; ++i)
    {
      for (const auto& pair : bpp[i])
      {
        uint j = pair.first;
        float prob = pair.second;
        
        // Only store probabilities above threshold
        if (prob > threshold() && j > i)
        {
          bp[i].push_back(std::make_pair(j, prob));
        }
      }
    }
  } catch (const std::exception& e) {
    // If LinFold fails, return empty BP
    bp.clear();
    bp.resize(seq.size());
  }
}

void 
LinFoldWrapper::calculate(const std::string& seq, const std::string& str, BP& bp)
{
  try {
    std::vector<std::vector<std::pair<u_int32_t, float>>> bpp;
    
    // Convert structure constraint to LinFold format
    // In DAFS: '.' = unpaired, '(' ')' = paired, '?' = any
    // In LinFold constraint format, we need to set up the constraint structure
    std::vector<uint32_t> constraint(str.size() + 1, _Fold::Options::ANY);
    
    // Parse the structure string to build constraints
    std::vector<int> stack;
    for (size_t i = 0; i < str.size(); ++i)
    {
      if (str[i] == '(')
      {
        stack.push_back(i);
      }
      else if (str[i] == ')')
      {
        if (!stack.empty())
        {
          int j = stack.back();
          stack.pop_back();
          // LinFold uses 1-based indexing for constraints
          constraint[j + 1] = i + 1;
          constraint[i + 1] = j + 1;
        }
      }
      else if (str[i] == '.')
      {
        constraint[i + 1] = _Fold::Options::UNPAIRED;
      }
      // '?' remains as ANY
    }
    
    if (model_type_ == ModelType::LPV) {
      // Initialize LinFold with TurnerNearestNeighbor parameters for this sequence
      auto turner_params = std::make_unique<TurnerNearestNeighbor>(seq);
      linfold_turner_ = std::make_unique<LinFold<TurnerNearestNeighbor>>(std::move(turner_params));
      
      // Set up options with beam size and constraint
      LinFold<TurnerNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      opt.constraints(constraint);
      
      // First compute inside algorithm (forward pass)
      linfold_turner_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_turner_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities with constraints
      bpp = linfold_turner_->compute_basepairing_probabilities(seq, opt);
    } else { // ModelType::LPC
      // Initialize LinFold with CONTRAfoldNearestNeighbor parameters for this sequence
      auto contra_params = std::make_unique<CONTRAfoldNearestNeighbor>(seq);
      linfold_contra_ = std::make_unique<LinFold<CONTRAfoldNearestNeighbor>>(std::move(contra_params));
      
      // Set up options with beam size and constraint
      LinFold<CONTRAfoldNearestNeighbor>::Options opt;
      opt.beam_size(beam_size_);
      opt.constraints(constraint);
      
      // First compute inside algorithm (forward pass)
      linfold_contra_->compute_inside(seq, opt);
      
      // Then compute outside algorithm (backward pass)
      linfold_contra_->compute_outside(seq, opt);
      
      // Now compute base-pairing probabilities with constraints
      bpp = linfold_contra_->compute_basepairing_probabilities(seq, opt);
    }
    
    // Convert LinFold output format to BP format used by DAFS
    uint L = seq.size();
    bp.clear();
    bp.resize(L);
    
    for (uint i = 0; i < L; ++i)
    {
      for (const auto& pair : bpp[i])
      {
        uint j = pair.first;
        float prob = pair.second;
        
        // Only store probabilities above threshold
        if (prob > threshold() && j > i)
        {
          bp[i].push_back(std::make_pair(j, prob));
        }
      }
    }
  } catch (const std::exception& e) {
    // If LinFold fails, return empty BP
    bp.clear();
    bp.resize(seq.size());
  }
}