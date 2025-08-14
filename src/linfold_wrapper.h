/*
 * LinFold wrapper for DAFS
 * 
 * This file provides a wrapper class for LinFold that implements
 * the Fold::Model interface, allowing LinFold to be used
 * consistently with other folding algorithms in DAFS.
 */

#ifndef __INC_LINFOLD_WRAPPER_H__
#define __INC_LINFOLD_WRAPPER_H__

#include "fold.h"
#include "linearfold/param/turner.h"
#include "linearfold/param/contrafold.h"
#include "linearfold/fold/linfold.h"
#include <memory>

class LinFoldWrapper : public Fold::Model
{
public:
  enum class ModelType {
    LPV,  // LinearPartition-V (Turner model)
    LPC   // LinearPartition-C (CONTRAfold model)
  };
  
  LinFoldWrapper(float th, ModelType model_type = ModelType::LPV, uint32_t beam_size = 100);
  virtual ~LinFoldWrapper() = default;
  
  virtual void calculate(const std::string& seq, BP& bp) override;
  virtual void calculate(const std::string& seq, const std::string& str, BP& bp) override;
  
  void set_beam_size(uint32_t beam_size) { beam_size_ = beam_size; }
  uint32_t get_beam_size() const { return beam_size_; }
  ModelType get_model_type() const { return model_type_; }

private:
  uint32_t beam_size_;
  ModelType model_type_;
  std::unique_ptr<LinFold<TurnerNearestNeighbor>> linfold_turner_;
  std::unique_ptr<LinFold<CONTRAfoldNearestNeighbor>> linfold_contra_;
};

#endif // __INC_LINFOLD_WRAPPER_H__