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

GradientManager::GradientManager(Method method, float eta0, float lb)
    : method_(method), eta0_(eta0), lb_(lb), current_eta_(eta0), 
      violations_(0), step_count_(0.0)
{
}

void GradientManager::initialize(uint L1, uint L2)
{
    // Initialize Lagrange multipliers
    q_x_.assign(L1, VF(L1, 0.0));
    q_y_.assign(L2, VF(L2, 0.0));
    q_z_.assign(L1, VF(L2, 0.0));
    
    // Initialize gradient history for adaptive methods
    if (method_ == ADAGRAD) {
        g2_x_.assign(L1, VF(L1, 0.0));
        g2_y_.assign(L2, VF(L2, 0.0));
        g2_z_.assign(L1, VF(L2, 0.0));
    }
    else if (method_ == ADAM) {
        m_x_.assign(L1, VF(L1, 0.0));
        v_x_.assign(L1, VF(L1, 0.0));
        m_y_.assign(L2, VF(L2, 0.0));
        v_y_.assign(L2, VF(L2, 0.0));
        m_z_.assign(L1, VF(L2, 0.0));
        v_z_.assign(L1, VF(L2, 0.0));
    }
}

void GradientManager::set_multipliers(const VVF& q_x, const VVF& q_y, const VVF& q_z)
{
    q_x_ = q_x;
    q_y_ = q_y;
    q_z_ = q_z;
}

void GradientManager::get_multipliers(VVF& q_x, VVF& q_y, VVF& q_z) const
{
    q_x = q_x_;
    q_y = q_y_;
    q_z = q_z_;
}

float GradientManager::adagrad_update(float& g2, float grad)
{
    const float eps = 1e-6;
    g2 += grad * grad;
    return eta0_ * grad / std::sqrt(g2 + eps);
}

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

void GradientManager::update_standard_stepsize(float score, float prev_score, uint cbp_size, uint t)
{
    if (score > prev_score || t == 0) {
        if (cbp_size > 0) {
            step_count_ += std::max(0.0f, 4.0f * cbp_size - violations_) / (4.0f * cbp_size);
            current_eta_ = eta0_ / (1.0 + step_count_);
        }
    }
}

uint GradientManager::update_gradients(const std::vector<CBP>& cbp,
                                      const VU& x, const VU& y, const VU& z,
                                      const VVU& c_x, const VVU& c_y, const VVU& c_z,
                                      uint t, float score, float prev_score)
{
    const uint L1 = q_x_.size();
    const uint L2 = q_y_.size();
    
    // Reset violation count
    violations_ = 0;
    
    // Compute constraint violations and gradients
    VVI t_x(L1, VI(L1, 0));
    VVI t_y(L2, VI(L2, 0));
    VVI t_z(L1, VI(L2, 0));
    
    // Check consensus base-pair constraints
    for (uint u = 0; u != cbp.size(); ++u) {
        const uint i = cbp[u].first.first, j = cbp[u].first.second;
        const uint k = cbp[u].second.first, l = cbp[u].second.second;
        const float s_w = q_x_[i][j] + q_y_[k][l] - q_z_[i][k] - q_z_[j][l];
        
        if (s_w > 0.0f) {
            t_x[i][j]++;
            t_y[k][l]++;
            t_z[i][k]++;
            t_z[j][l]++;
        }
    }
    
    // Update Lagrangian for x (=q_x) using sparse update
    for (uint i = 0; i != L1; ++i) {
        const uint j = x[i];
        if (j != -1u && t_x[i][j] != 1) {
            violations_++;
            float grad = t_x[i][j] - 1;
            
            switch (method_) {
                case ADAGRAD:
                    q_x_[i][j] -= adagrad_update(g2_x_[i][j], grad);
                    break;
                case ADAM:
                    q_x_[i][j] -= adam_update(m_x_[i][j], v_x_[i][j], grad, t + 1);
                    break;
                case STANDARD:
                    q_x_[i][j] -= current_eta_ * grad;
                    break;
            }
        }
        
        for (uint jj = 0; jj != c_x[i].size(); ++jj) {
            const uint j = c_x[i][jj];
            if (x[i] != j && t_x[i][j] != 0) {
                violations_++;
                float grad = t_x[i][j];
                
                switch (method_) {
                    case ADAGRAD:
                        q_x_[i][j] -= adagrad_update(g2_x_[i][j], grad);
                        break;
                    case ADAM:
                        q_x_[i][j] -= adam_update(m_x_[i][j], v_x_[i][j], grad, t + 1);
                        break;
                    case STANDARD:
                        q_x_[i][j] -= current_eta_ * grad;
                        break;
                }
            }
        }
    }
    
    // Update Lagrangian for y (=q_y) using sparse update
    for (uint k = 0; k != L2; ++k) {
        const uint l = y[k];
        if (l != -1u && t_y[k][l] != 1) {
            violations_++;
            float grad = t_y[k][l] - 1;
            
            switch (method_) {
                case ADAGRAD:
                    q_y_[k][l] -= adagrad_update(g2_y_[k][l], grad);
                    break;
                case ADAM:
                    q_y_[k][l] -= adam_update(m_y_[k][l], v_y_[k][l], grad, t + 1);
                    break;
                case STANDARD:
                    q_y_[k][l] -= current_eta_ * grad;
                    break;
            }
        }
        
        for (uint ll = 0; ll != c_y[k].size(); ++ll) {
            const uint l = c_y[k][ll];
            if (y[k] != l && t_y[k][l] != 0) {
                violations_++;
                float grad = t_y[k][l];
                
                switch (method_) {
                    case ADAGRAD:
                        q_y_[k][l] -= adagrad_update(g2_y_[k][l], grad);
                        break;
                    case ADAM:
                        q_y_[k][l] -= adam_update(m_y_[k][l], v_y_[k][l], grad, t + 1);
                        break;
                    case STANDARD:
                        q_y_[k][l] -= current_eta_ * grad;
                        break;
                }
            }
        }
    }
    
    // Update Lagrangian for z (=q_z) using sparse update
    for (uint i = 0; i != L1; ++i) {
        const uint k = z[i];
        if (k != -1u) {
            if (t_z[i][k] > 1) {
                violations_++;
            }
            float grad = 1 - t_z[i][k];
            float update = 0.0;
            
            switch (method_) {
                case ADAGRAD:
                    update = adagrad_update(g2_z_[i][k], grad);
                    break;
                case ADAM:
                    update = adam_update(m_z_[i][k], v_z_[i][k], grad, t + 1);
                    break;
                case STANDARD:
                    update = current_eta_ * grad;
                    break;
            }
            
            q_z_[i][k] = std::max(0.0f, q_z_[i][k] - update);
        }
        
        for (uint kk = 0; kk != c_z[i].size(); ++kk) {
            const uint k = c_z[i][kk];
            if (z[i] != k) {
                if (t_z[i][k] > 0) {
                    violations_++;
                }
                float grad = -t_z[i][k];
                float update = 0.0;
                
                switch (method_) {
                    case ADAGRAD:
                        update = adagrad_update(g2_z_[i][k], grad);
                        break;
                    case ADAM:
                        update = adam_update(m_z_[i][k], v_z_[i][k], grad, t + 1);
                        break;
                    case STANDARD:
                        update = current_eta_ * grad;
                        break;
                }
                
                q_z_[i][k] = std::max(0.0f, q_z_[i][k] - update);
            }
        }
    }
    
    // Update step size for standard method
    if (method_ == STANDARD) {
        update_standard_stepsize(score, prev_score, cbp.size(), t);
    }
    
    return violations_;
}