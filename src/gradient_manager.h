/*
 * Copyright (C) 2024 DAFS Developers
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

#ifndef __INC_GRADIENT_MANAGER_H__
#define __INC_GRADIENT_MANAGER_H__

#include "typedefs.h"
#include "sparse_matrix.h"
#include <map>
#include <memory>
#include <cmath>

//#define USE_ADAGRAD
//#define USE_ADAM
#define USE_ADAPTIVE

// Consensus base-pair indices
typedef std::pair<std::pair<uint, uint>, std::pair<uint, uint>> CBP;


class GradientManager {
public:
    GradientManager(float eta0, float lb = 0.0, float gradient_clip = 1.0,
                    bool sparse_structure_storage = false,
                    bool sparse_alignment_storage = false);
    ~GradientManager() = default;
    
    // Initialize gradient matrices
    void initialize(uint L1, uint L2);
    
    // Set/Get Lagrange multipliers
    void set_multipliers(const VVF& q_x, const VVF& q_y, const VVF& q_z);
    void get_multipliers(VVF& q_x, VVF& q_y, VVF& q_z) const;
    bool uses_sparse_structure_storage() const { return sparse_structure_storage_; }
    bool uses_sparse_alignment_storage() const { return sparse_alignment_storage_; }
    const VVF& dense_q_x() const { return q_x_; }
    const VVF& dense_q_y() const { return q_y_; }
    const VVF& dense_q_z() const { return q_z_; }
    const SparseFloatMatrix& sparse_q_x() const { return sparse_q_x_; }
    const SparseFloatMatrix& sparse_q_y() const { return sparse_q_y_; }
    const SparseFloatMatrix& sparse_q_z() const { return sparse_q_z_; }
    float q_x(uint i, uint j) const;
    float q_y(uint i, uint j) const;
    float q_z(uint i, uint j) const;
    
    // Update gradients based on constraint violations
    uint update_gradients(const std::vector<CBP>& cbp,
                         const VU& x, const VU& y, const VU& z, const VU& w_cbp,
                         const VVU& c_x, const VVU& c_y, const VVU& c_z,
                         uint t, float score);
    
    // Get current violation count
    uint get_violations() const { return violations_; }
    
    // Get the current Polyak scaling coefficient alpha_t.  The actual step is
    // alpha_t * (L(q_t) - LB) / ||g_t||^2.
    float get_step_size() const { return current_eta_; }
    // Actual multiplier applied to the materialized subgradient in the most
    // recent update.  Unlike get_step_size(), this includes the duality gap
    // and squared-gradient denominator.
    float get_last_update_size() const { return last_update_size_; }
    
    // Set lower bound for adaptive method
    void set_lower_bound(float lb) { lb_ = lb; }
    
    // Set gradient clipping threshold
    void set_gradient_clip(float clip) { gradient_clip_ = clip; }
    
    
private:
    float eta0_;          // Initial step size
    float lb_;            // Lower bound
    float current_eta_;   // Current step size
    float last_update_size_; // Actual Polyak multiplier used most recently
    uint violations_;     // Number of constraint violations
    float gradient_clip_; // Gradient clipping threshold
    bool sparse_structure_storage_; // Sparse q_x/q_y for LinearFold models
    bool sparse_alignment_storage_; // Sparse q_z for LinearAlign
    
    // Lagrange multipliers
    VVF q_x_, q_y_, q_z_;
    SparseFloatMatrix sparse_q_x_, sparse_q_y_, sparse_q_z_;
    
    // For adaptive methods
#if defined(USE_ADAGRAD)
    VVF g2_x_, g2_y_, g2_z_;                    // AdaGrad: accumulated squared gradients
    SparseFloatMatrix sparse_g2_x_, sparse_g2_y_, sparse_g2_z_;
#elif defined(USE_ADAM)
    VVF m_x_, m_y_, v_x_, v_y_, m_z_, v_z_;    // Adam: first and second moment estimates
    SparseFloatMatrix sparse_m_x_, sparse_m_y_, sparse_v_x_, sparse_v_y_;
    SparseFloatMatrix sparse_m_z_, sparse_v_z_;
#endif
    
    // Helper functions for different update methods
#if defined(USE_ADAGRAD)
    float adagrad_update(float& g2, float grad);
#elif defined(USE_ADAM)
    float adam_update(float& m, float& v, float grad, uint t);
#endif
    void update_standard_stepsize(float score, float prev_score, uint cbp_size, uint t);
    void update_adaptive_stepsize(float score, const VU& x, const VU& y, const VU& z,
                                 const VVU& c_x, const VVU& c_y, const VVU& c_z);
    void compute_constraint_violations(const std::vector<CBP>& cbp);
    
    // Gradient clipping helper
    float clip_update(float update) const;
};

#endif // __INC_GRADIENT_MANAGER_H__
