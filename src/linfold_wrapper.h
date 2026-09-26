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
#include "linearfold/param/profile.h"
#include "linearfold/fold/linfold.h"
#include <memory>
#include <mutex>

class LinFoldWrapper : public Fold::Model
{
public:
  enum class ModelType {
    LPV,  // LinearPartition-V (Turner model)
    LPC   // LinearPartition-C (CONTRAfold model)
  };

  enum class ProfileEnergyMode {
    Legacy,
    RNAalifold
  };
  
  LinFoldWrapper(float th, ModelType model_type = ModelType::LPV,
                 uint32_t beam_size = 100,
                 ProfileEnergyMode profile_energy_mode = ProfileEnergyMode::Legacy);
  virtual ~LinFoldWrapper() = default;
  
  virtual void calculate(const std::string& seq, BP& bp) override;
  virtual void calculate(const std::string& seq, const std::string& str, BP& bp) override;
  // Consensus partition function and BPPs for an aligned profile.  The
  // optional constraint uses DAFS's '.', '?', '(', ')' alphabet.
  // Returns log(Z) at 37 C; profile scores use the thermodynamic model for
  // both LPV and LPC, just as RNAalifold is independent of the row predictor.
  // Internal DAFS refinement may mark unsupported predicted pairs unpaired;
  // explicit caller constraints remain strict by default.
  double calculate_profile(const ALN& alignment, const std::vector<Fasta>& sequences,
                           BP& bp, const std::string& constraint = "",
                           bool relax_unsupported_pairs = false) const;
  
  void set_beam_size(uint32_t beam_size) { beam_size_ = beam_size; }
  uint32_t get_beam_size() const { return beam_size_; }
  ModelType get_model_type() const { return model_type_; }
  ProfileEnergyMode get_profile_energy_mode() const { return profile_energy_mode_; }

private:
  uint32_t beam_size_;
  ModelType model_type_;
  ProfileEnergyMode profile_energy_mode_;
  std::unique_ptr<LinFold<TurnerNearestNeighbor>> linfold_turner_;
  std::unique_ptr<LinFold<CONTRAfoldNearestNeighbor>> linfold_contra_;
  mutable std::mutex profile_mutex_;
  mutable std::unique_ptr<LinFold<ProfileNearestNeighbor>> linfold_profile_;
};

#endif // __INC_LINFOLD_WRAPPER_H__
