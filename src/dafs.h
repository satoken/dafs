/*
 * Copyright (C) 2012 Kengo Sato
 *
 * This file is part of DAFS.
 *
 * DAFS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * DAFS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with DAFS.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef DAFS_H
#define DAFS_H

#include <vector>
#include <memory>
#include <iostream>
#include <random>
#include <unordered_set>
#include "typedefs.h"
#include "fa.h"
#include "fold.h"
#include "align.h"
#include "gradient_manager.h"


class DAFS
{
private:
  // nodes in the guide tree
  typedef std::pair<float, std::pair<uint, uint>> node_t;

  struct CBPHash {
    size_t operator()(const CBP& candidate) const noexcept {
      const auto& [ij, kl] = candidate;
      size_t seed = std::hash<uint>{}(ij.first);
      seed ^= std::hash<uint>{}(ij.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= std::hash<uint>{}(kl.first) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= std::hash<uint>{}(kl.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };

public:
  DAFS();
  ~DAFS();

  DAFS &parse_options(int &argc, char **&argv);
  int run();

private:
  void relax_matching_probability();
  void relax_basepairing_probability();
  void relax_fourway_consistency();
  void build_tree();
  void print_tree(std::ostream &os, int i) const;
  void project_alignment(ALN &aln, const ALN &aln1, const ALN &aln2, const VU &z) const;
  void project_secondary_structure(VU &xx, VU &yy, const VU &x, const VU &y, const VU &z) const;
  void average_matching_probability(SparseFloatMatrix &posterior, const ALN &aln1, const ALN &aln2) const;
  void average_basepairing_probability(SparseFloatMatrix &posterior, const ALN &aln, bool use_alifold) const;
  void update_basepairing_probability(SparseFloatMatrix &posterior, const VU &ss, const std::string &str,
                                      const ALN &aln, bool use_alifold) const;
  std::string profile_consensus_sequence(const ALN &aln) const;
  void calculate_profile_basepairing_probability(const ALN &aln, BP &bp) const;
  void calculate_profile_basepairing_probability(const ALN &aln,
                                                  const std::string &constraint,
                                                  BP &bp) const;
  void align_alignments(ALN &aln, const ALN &aln1, const ALN &aln2);
  float align_alignments(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2);
  float calculate_alignment_only_score(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2) const;
  float repair_feasible_solution(VU& repaired_x, VU& repaired_y, VU& repaired_z,
                                 const VU& x, const VU& y, const VU& z,
                                 const SparseFloatMatrix& p_x,
                                 const SparseFloatMatrix& p_y,
                                 const SparseFloatMatrix& p_z,
                                 uint N1, uint N2, float min_th_s) const;
  float solve(VU &x, VU &y, VU &z, const SparseFloatMatrix &p_x,
              const SparseFloatMatrix &p_y, const SparseFloatMatrix &p_z,
              const ALN &aln1, const ALN &aln2);
  float solve_by_dd(VU &x, VU &y, VU &z, const SparseFloatMatrix &p_x,
                    const SparseFloatMatrix &p_y, const SparseFloatMatrix &p_z,
                    const ALN &aln1, const ALN &aln2);
  float solve_by_ip(VU &x, VU &y, VU &z, const SparseFloatMatrix &p_x,
                    const SparseFloatMatrix &p_y, const SparseFloatMatrix &p_z,
                    const ALN &aln1, const ALN &aln2) const;
  void align(ALN &aln, int ch);
  float align(VU &ss, ALN &aln, int ch);
  float refine(VU &ss, ALN &aln);
  void output_verbose(const VU &x, const VU &y, const VU &z, const ALN &aln1, const ALN &aln2) const;
  void output(std::ostream &os, const ALN &aln) const;
  void output(std::ostream &os, ALN::const_iterator b, ALN::const_iterator e) const;
  
  // Dynamic CBP generation methods
  bool is_valid_cbp(uint i, uint j, uint k, uint l, 
                    const SparseFloatMatrix& p_x,
                    const SparseFloatMatrix& p_y,
                    const SparseFloatMatrix& p_z,
                    uint N1, uint N2, float min_th_s) const;
  void add_cbp_if_new(const CBP& candidate, std::vector<CBP>& cbp, 
                      VVU& c_x, VVU& c_y, VVU& c_z);
  void generate_cbp_from_solution(const VU& x, const VU& y, const VU& z,
                                  const SparseFloatMatrix& p_x,
                                  const SparseFloatMatrix& p_y,
                                  const SparseFloatMatrix& p_z,
                                  uint N1, uint N2, float min_th_s,
                                  std::vector<CBP>& cbp, VVU& c_x, VVU& c_y, VVU& c_z);
  void generate_positive_reduced_cost_cbp(const GradientManager& gm,
                                          const SparseFloatMatrix& p_x,
                                          const SparseFloatMatrix& p_y,
                                          const SparseFloatMatrix& p_z,
                                          const VVU& p_z_forward,
                                          const VVU& p_z_reverse,
                                          const std::vector<std::pair<uint, uint>>& q_x_support,
                                          const std::vector<std::pair<uint, uint>>& q_y_support,
                                          uint N1, uint N2, float min_th_s,
                                          std::vector<CBP>& cbp,
                                          VVU& c_x, VVU& c_y, VVU& c_z);

private:
  float w_pct_a_;                   // the weight of PCT for alignment matching probabilities
  float w_pct_s_;                   // the weight of PCT for base-pairing probabilities
  float w_pct_f_;                   // the weight of four-way PCT
  uint n_refinement_;               // the number of the iterative refinement
  uint t_max_;                      // the maximum number of the iteration of the subgradient update
  float th_a_;                      // the threshold for base-pairing probabilities
  VF th_s_;                         // the threshold for alignment matching probabilities
  float w_;                         // the weight for base pairs in the objective function
  float eta0_;                      // the initial step width of the subgradient update
  std::unique_ptr<Align::Model> a_model_;           // alignment model
  std::unique_ptr<Align::Decoder> a_decoder_;       // alignment decoder
  std::unique_ptr<Fold::Model> s_model_;            // folding model
  std::unique_ptr<Fold::Decoder> s_decoder_;        // folding decoder
  std::unique_ptr<Fold::Decoder> s_decoder1_;       // folding decoder for the final folding
  std::vector<Fasta> fa_;           // input sequences
  std::vector<std::vector<MP>> mp_; // alignment matching probability matrices
  std::vector<BP> bp_;              // base-pairing probability matrices
  VVF sim_;                         // simalarity matrix between input sequences
  std::vector<node_t> tree_;        // guide tree
  
  // Dynamic CBP generation
  bool use_dynamic_cbp_;            // whether to use dynamic CBP generation
  bool use_sparse_structure_lagrangian_; // sparse q_x/q_y for LinearFold
  bool use_sparse_alignment_lagrangian_; // sparse q_z for LinearAlign
  bool use_linear_structure_decoder_;    // beam max decoder for folding
  bool use_linear_alignment_decoder_;    // beam max decoder for alignment
  std::unordered_set<CBP, CBPHash> cbp_set_; // for expected O(1) duplicate removal
  
  bool use_alifold_;
  bool use_alifold1_;
  bool use_linear_profile_folding_;
  bool use_bp_update_;
  bool use_bp_update1_;
  // bool use_bpscore_;
  uint verbose_;
  mutable std::mt19937 g_; // random number generator for shuffling
};

#endif // DAFS_H
