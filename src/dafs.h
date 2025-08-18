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
#include "typedefs.h"
#include "fa.h"
#include "fold.h"
#include "align.h"


class DAFS
{
private:
  // nodes in the guide tree
  typedef std::pair<float, std::pair<uint, uint>> node_t;

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
  void average_matching_probability(VVF &posterior, const ALN &aln1, const ALN &aln2) const;
  void average_basepairing_probability(VVF &posterior, const ALN &aln, bool use_alifold) const;
  void update_basepairing_probability(VVF &posterior, const VU &ss, const std::string &str,
                                      const ALN &aln, bool use_alifold) const;
  void align_alignments(ALN &aln, const ALN &aln1, const ALN &aln2);
  float align_alignments(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2);
  float calculate_alignment_only_score(VU &ss, ALN &aln, const ALN &aln1, const ALN &aln2) const;
  float solve(VU &x, VU &y, VU &z, const VVF &p_x, const VVF &p_y, const VVF &p_z,
              const ALN &aln1, const ALN &aln2);
  float solve_by_dd(VU &x, VU &y, VU &z, const VVF &p_x, const VVF &p_y, const VVF &p_z,
                    const ALN &aln1, const ALN &aln2);
  float solve_by_ip(VU &x, VU &y, VU &z, const VVF &p_x, const VVF &p_y, const VVF &p_z,
                    const ALN &aln1, const ALN &aln2) const;
  void align(ALN &aln, int ch);
  float align(VU &ss, ALN &aln, int ch);
  float refine(VU &ss, ALN &aln);
  void output_verbose(const VU &x, const VU &y, const VU &z, const ALN &aln1, const ALN &aln2) const;
  void output(std::ostream &os, const ALN &aln) const;
  void output(std::ostream &os, ALN::const_iterator b, ALN::const_iterator e) const;

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
  bool use_alifold_;
  bool use_alifold1_;
  bool use_bp_update_;
  bool use_bp_update1_;
  // bool use_bpscore_;
  uint verbose_;
  mutable std::mt19937 g_; // random number generator for shuffling
};

#endif // DAFS_H