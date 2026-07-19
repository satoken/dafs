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

#include "gradient_manager.h"
#include <algorithm>
#include <cassert>
#include <vector>

#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/stopwatch.h"

namespace {

enum class GradientBlock {
    X,
    Y,
    Z,
};

struct GradientEntry {
    GradientBlock block;
    uint row;
    uint column;
    float value;
    bool violated;
};

} // namespace

GradientManager::GradientManager(float eta0, float lb, float gradient_clip,
                                 bool sparse_structure_storage,
                                 bool sparse_alignment_storage)
    : eta0_(eta0), lb_(lb), current_eta_(eta0), 
      violations_(0), gradient_clip_(gradient_clip),
      sparse_structure_storage_(sparse_structure_storage),
      sparse_alignment_storage_(sparse_alignment_storage)
{
}

void GradientManager::initialize(uint L1, uint L2)
{
    // Initialize Lagrange multipliers in exactly one representation.
    if (sparse_structure_storage_) {
        sparse_q_x_.assign(L1, L1);
        sparse_q_y_.assign(L2, L2);
    } else {
        q_x_.assign(L1, VF(L1, 0.0));
        q_y_.assign(L2, VF(L2, 0.0));
    }
    if (sparse_alignment_storage_) {
        sparse_q_z_.assign(L1, L2);
    } else {
        q_z_.assign(L1, VF(L2, 0.0));
    }
    
    // Initialize gradient history for adaptive methods
#if defined(USE_ADAGRAD)
    if (sparse_structure_storage_) {
        sparse_g2_x_.assign(L1, L1);
        sparse_g2_y_.assign(L2, L2);
    } else {
        g2_x_.assign(L1, VF(L1, 0.0));
        g2_y_.assign(L2, VF(L2, 0.0));
    }
    if (sparse_alignment_storage_)
        sparse_g2_z_.assign(L1, L2);
    else
        g2_z_.assign(L1, VF(L2, 0.0));
#elif defined(USE_ADAM)
    if (sparse_structure_storage_) {
        sparse_m_x_.assign(L1, L1); sparse_v_x_.assign(L1, L1);
        sparse_m_y_.assign(L2, L2); sparse_v_y_.assign(L2, L2);
    } else {
        m_x_.assign(L1, VF(L1, 0.0)); v_x_.assign(L1, VF(L1, 0.0));
        m_y_.assign(L2, VF(L2, 0.0)); v_y_.assign(L2, VF(L2, 0.0));
    }
    if (sparse_alignment_storage_) {
        sparse_m_z_.assign(L1, L2); sparse_v_z_.assign(L1, L2);
    } else {
        m_z_.assign(L1, VF(L2, 0.0)); v_z_.assign(L1, VF(L2, 0.0));
    }
#endif
}

void GradientManager::set_multipliers(const VVF& q_x, const VVF& q_y, const VVF& q_z)
{
    assert(!sparse_structure_storage_ && !sparse_alignment_storage_);
    q_x_ = q_x;
    q_y_ = q_y;
    q_z_ = q_z;
}

void GradientManager::get_multipliers(VVF& q_x, VVF& q_y, VVF& q_z) const
{
    assert(!sparse_structure_storage_ && !sparse_alignment_storage_);
    q_x = q_x_;
    q_y = q_y_;
    q_z = q_z_;
}

float GradientManager::q_x(uint i, uint j) const
{
    return sparse_structure_storage_ ? sparse_q_x_.get(i, j) : q_x_[i][j];
}

float GradientManager::q_y(uint i, uint j) const
{
    return sparse_structure_storage_ ? sparse_q_y_.get(i, j) : q_y_[i][j];
}

float GradientManager::q_z(uint i, uint j) const
{
    return sparse_alignment_storage_ ? sparse_q_z_.get(i, j) : q_z_[i][j];
}

#if defined(USE_ADAGRAD)
float GradientManager::adagrad_update(float& g2, float grad)
{
    const float eps = 1e-6;
    g2 += grad * grad;
    return eta0_ * grad / std::sqrt(g2 + eps);
}
#endif

#if defined(USE_ADAM)
float GradientManager::adam_update(float& m, float& v, float grad, uint t)
{
    const float beta1 = 0.9;
    const float beta2 = 0.999;
    const float eps = 1e-8;
    
    m = beta1 * m + (1 - beta1) * grad;
    v = beta2 * v + (1 - beta2) * grad * grad;
    
    const float m_hat = m / (1 - std::pow(beta1, t));
    const float v_hat = v / (1 - std::pow(beta2, t));
    
    return eta0_ * m_hat / (std::sqrt(v_hat) + eps);
}
#endif

uint GradientManager::update_gradients(const std::vector<CBP>& cbp,
                                      const VU& x, const VU& y, const VU& z, const VU& w_cbp,
                                      const VVU& c_x, const VVU& c_y, const VVU& c_z,
                                      uint t, float score)
{
    const uint L1 = sparse_structure_storage_ ? sparse_q_x_.rows() : q_x_.size();
    const uint L2 = sparse_structure_storage_ ? sparse_q_y_.rows() : q_y_.size();
    
    // Reset violation count
    violations_ = 0;
    
    assert(c_x.size() == L1 && c_y.size() == L2 && c_z.size() == L1);

    // Store consensus counts in arrays aligned with the sorted projection
    // lists.  Their total size is |c_x|+|c_y|+|c_z| instead of the three
    // dense Cartesian products.
    VVI t_x(L1), t_y(L2), t_z(L1);
    for (uint i = 0; i != L1; ++i) {
        t_x[i].assign(c_x[i].size(), 0);
        t_z[i].assign(c_z[i].size(), 0);
    }
    for (uint k = 0; k != L2; ++k)
        t_y[k].assign(c_y[k].size(), 0);

    const auto projection_index = [](const VU& projection, uint column) {
        const auto it = std::lower_bound(projection.begin(), projection.end(), column);
        assert(it != projection.end() && *it == column);
        return static_cast<size_t>(it - projection.begin());
    };

    const auto projection_count = [](const VU& projection, const VI& counts,
                                     uint column) {
        const auto it = std::lower_bound(projection.begin(), projection.end(), column);
        return it == projection.end() || *it != column
             ? 0
             : counts[static_cast<size_t>(it - projection.begin())];
    };
    
    // Check consensus base-pair constraints
    for (const auto& u: w_cbp) {
        // w_ijkl=1
        const auto &[i, j] = cbp[u].first;
        const auto &[k, l] = cbp[u].second;
        ++t_x[i][projection_index(c_x[i], j)];
        ++t_y[k][projection_index(c_y[k], l)];
        ++t_z[i][projection_index(c_z[i], k)];
        ++t_z[j][projection_index(c_z[j], l)];
    }
   
    // Materialize the sparse subgradient once.  The exact same coordinates
    // are used below for the norm, violation count, and multiplier update.
    // This prevents the Polyak denominator from depending on a different
    // compile-time update configuration.
    std::vector<GradientEntry> gradients;
    gradients.reserve(w_cbp.size() * 4 + L1 + L2);

    for (uint i = 0; i != L1; ++i) {
        const uint j = x[i];
        const int selected_count = j == -1u
                                 ? 0
                                 : projection_count(c_x[i], t_x[i], j);
        if (j != -1u && selected_count != 1) {
            gradients.push_back(
                {GradientBlock::X, i, j,
                 static_cast<float>(selected_count - 1), true}); // x_ij=1
        }
        
        for (size_t v = 0; v != c_x[i].size(); ++v) {
            const uint candidate = c_x[i][v];
            if (x[i] != candidate && t_x[i][v] != 0) {
                gradients.push_back(
                    {GradientBlock::X, i, candidate,
                     static_cast<float>(t_x[i][v]), true}); // x_ij=0
            }
        }
    }
    
    for (uint k = 0; k != L2; ++k) {
        const uint l = y[k];
        const int selected_count = l == -1u
                                 ? 0
                                 : projection_count(c_y[k], t_y[k], l);
        if (l != -1u && selected_count != 1) {
            gradients.push_back(
                {GradientBlock::Y, k, l,
                 static_cast<float>(selected_count - 1), true}); // y_kl=1
        }
        
        for (size_t v = 0; v != c_y[k].size(); ++v) {
            const uint candidate = c_y[k][v];
            if (y[k] != candidate && t_y[k][v] != 0) {
                gradients.push_back(
                    {GradientBlock::Y, k, candidate,
                     static_cast<float>(t_y[k][v]), true}); // y_kl=0
            }
        }
    }
    
    for (uint i = 0; i != L1; ++i) {
        const uint k = z[i];
        if (k != -1u) {
            const int selected_count = projection_count(c_z[i], t_z[i], k);
            const float grad = 1 - selected_count; // z_ik=1
            gradients.push_back(
                {GradientBlock::Z, i, k, grad, grad < 0.0f});
        }
        
        for (size_t v = 0; v != c_z[i].size(); ++v) {
            const uint candidate = c_z[i][v];
            if (z[i] != candidate) {
                const float grad = -t_z[i][v]; // z_ik=0
                gradients.push_back(
                    {GradientBlock::Z, i, candidate, grad, grad < 0.0f});
            }
        }
    }

    float g2 = 0.0f;
    for (const GradientEntry& gradient : gradients)
        g2 += gradient.value * gradient.value;

    // Polyak step for minimizing the Lagrangian dual upper bound:
    //   eta_t = alpha_t * (L(q_t) - LB) / ||g_t||^2.
    // Clamp the duality gap to zero so rounding error cannot reverse the
    // subgradient direction.  With a zero subgradient no update is needed.
    const float duality_gap = std::max(0.0f, score - lb_);
    const float eta = g2 > 0.0f
                    ? current_eta_ * duality_gap / g2
                    : 0.0f;
    spdlog::debug("eta: {}, gap: {}, g^2: {}, score: {}, lb_: {}",
                  eta, duality_gap, g2, score, lb_);

    for (const GradientEntry& gradient : gradients) {
      if (gradient.violated)
          ++violations_;

      const bool use_sparse = gradient.block == GradientBlock::Z
                            ? sparse_alignment_storage_
                            : sparse_structure_storage_;
      if (use_sparse) {
        SparseFloatMatrix* multiplier;
#if defined(USE_ADAGRAD)
        SparseFloatMatrix* accumulated_g2;
#elif defined(USE_ADAM)
        SparseFloatMatrix* first_moment;
        SparseFloatMatrix* second_moment;
#endif

        switch (gradient.block) {
        case GradientBlock::X:
            multiplier = &sparse_q_x_;
#if defined(USE_ADAGRAD)
            accumulated_g2 = &sparse_g2_x_;
#elif defined(USE_ADAM)
            first_moment = &sparse_m_x_; second_moment = &sparse_v_x_;
#endif
            break;
        case GradientBlock::Y:
            multiplier = &sparse_q_y_;
#if defined(USE_ADAGRAD)
            accumulated_g2 = &sparse_g2_y_;
#elif defined(USE_ADAM)
            first_moment = &sparse_m_y_; second_moment = &sparse_v_y_;
#endif
            break;
        case GradientBlock::Z:
            multiplier = &sparse_q_z_;
#if defined(USE_ADAGRAD)
            accumulated_g2 = &sparse_g2_z_;
#elif defined(USE_ADAM)
            first_moment = &sparse_m_z_; second_moment = &sparse_v_z_;
#endif
            break;
        }

#if defined(USE_ADAGRAD)
        float history = accumulated_g2->get(gradient.row, gradient.column);
        const float update = clip_update(adagrad_update(history, gradient.value));
        accumulated_g2->set(gradient.row, gradient.column, history);
#elif defined(USE_ADAM)
        float m = first_moment->get(gradient.row, gradient.column);
        float v = second_moment->get(gradient.row, gradient.column);
        const float update = clip_update(adam_update(m, v, gradient.value, t + 1));
        first_moment->set(gradient.row, gradient.column, m);
        second_moment->set(gradient.row, gradient.column, v);
#else
        const float update = clip_update(eta * gradient.value);
#endif
        float value = multiplier->get(gradient.row, gradient.column) - update;
        if (gradient.block == GradientBlock::Z)
            value = std::max(0.0f, value);
        multiplier->set(gradient.row, gradient.column, value);
        continue;
      }

      float* multiplier;
#if defined(USE_ADAGRAD)
      float* accumulated_g2;
#elif defined(USE_ADAM)
      float* first_moment;
      float* second_moment;
#endif

        switch (gradient.block) {
        case GradientBlock::X:
            multiplier = &q_x_[gradient.row][gradient.column];
#if defined(USE_ADAGRAD)
            accumulated_g2 = &g2_x_[gradient.row][gradient.column];
#elif defined(USE_ADAM)
            first_moment = &m_x_[gradient.row][gradient.column];
            second_moment = &v_x_[gradient.row][gradient.column];
#endif
            break;
        case GradientBlock::Y:
            multiplier = &q_y_[gradient.row][gradient.column];
#if defined(USE_ADAGRAD)
            accumulated_g2 = &g2_y_[gradient.row][gradient.column];
#elif defined(USE_ADAM)
            first_moment = &m_y_[gradient.row][gradient.column];
            second_moment = &v_y_[gradient.row][gradient.column];
#endif
            break;
        case GradientBlock::Z:
            multiplier = &q_z_[gradient.row][gradient.column];
#if defined(USE_ADAGRAD)
            accumulated_g2 = &g2_z_[gradient.row][gradient.column];
#elif defined(USE_ADAM)
            first_moment = &m_z_[gradient.row][gradient.column];
            second_moment = &v_z_[gradient.row][gradient.column];
#endif
            break;
        }

#if defined(USE_ADAGRAD)
        const float update = clip_update(
            adagrad_update(*accumulated_g2, gradient.value));
#elif defined(USE_ADAM)
        const float update = clip_update(
            adam_update(*first_moment, *second_moment,
                        gradient.value, t + 1));
#else
        const float update = clip_update(eta * gradient.value);
#endif
        *multiplier -= update;
        if (gradient.block == GradientBlock::Z)
            *multiplier = std::max(0.0f, *multiplier);
    }
    
    return violations_;
}

float GradientManager::clip_update(float update) const {
    if (std::abs(update) > gradient_clip_) {
        return update > 0 ? gradient_clip_ : -gradient_clip_;
    }
    return update;
}
