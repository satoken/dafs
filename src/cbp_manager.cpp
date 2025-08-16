/*
 * CBP Manager implementation for Column Generation in DAFS
 */

#include "cbp_manager.h"
#include <algorithm>
#include <cassert>
#include <iostream>
#include "spdlog/spdlog.h"

CBPManager::CBPManager(uint L1, uint L2, float cutoff)
    : L1_(L1), L2_(L2), cutoff_(cutoff),
      p_x_(nullptr), p_y_(nullptr), p_z_(nullptr)
{
    c_x_.clear();
    c_x_.resize(L1_);
    c_y_.clear();
    c_y_.resize(L2_);
    c_z_.clear();
    c_z_.resize(L1_);
}

void CBPManager::initialize(const VVF& p_x, const VVF& p_y, const VVF& p_z,
                           uint N1, uint N2, float w, float th_a, 
                           const std::vector<float>& th_s)
{
    p_x_ = &p_x;
    p_y_ = &p_y;
    p_z_ = &p_z;
    N1_ = N1;
    N2_ = N2;
    w_ = w;
    th_a_ = th_a;
    th_s_ = th_s;
    
    // Clear existing data
    cbp_.clear();
    cbp_set_.clear();
    for (auto& c : c_x_) c.clear();
    for (auto& c : c_y_) c.clear();
    for (auto& c : c_z_) c.clear();
}

void CBPManager::generateInitialCBPs()
{
    if (!p_x_ || !p_y_ || !p_z_) {
        spdlog::error("CBPManager not initialized with probability matrices");
        return;
    }
    
    // True column generation: start with empty CBP set
    cbp_.clear();
    cbp_set_.clear();
    for (auto& c : c_x_) c.clear();
    for (auto& c : c_y_) c.clear();
    for (auto& c : c_z_) c.clear();
    
    std::cout << "CBPManager: Starting with empty CBP set (true column generation)" << std::endl;
    std::cout << "CBPManager: L1=" << L1_ << ", L2=" << L2_ << ", cutoff_=" << cutoff_ << ", w_=" << w_ << ", th_a_=" << th_a_ << std::endl;
    
    updateProjections();
}

bool CBPManager::shouldAddCBP(uint i, uint j, uint k, uint l) const
{
    if (!p_x_ || !p_y_ || !p_z_) return false;
    
    // Check basic probability cutoffs
    if ((*p_x_)[i][j] <= cutoff_ || (*p_y_)[k][l] <= cutoff_ ||
        (*p_z_)[i][k] <= cutoff_ || (*p_z_)[j][l] <= cutoff_) {
        return false;
    }
    
    // Check if already exists
    if (cbpExists(i, j, k, l)) {
        return false;
    }
    
    // For column generation, be extremely aggressive about adding CBPs to resolve violations
    // In true column generation, we must add CBPs to resolve violations, even if they seem weak
    float min_th_s = *std::min_element(th_s_.begin(), th_s_.end());
    float p = (N1_ * (*p_x_)[i][j] + N2_ * (*p_y_)[k][l]) / (N1_ + N2_);
    float q = ((*p_z_)[i][k] + (*p_z_)[j][l]) / 2;
    
    // Extremely aggressive - almost always accept if basic probabilities are above cutoff
    return (p > 0.0 && q > 0.0); // Only require positive probabilities
}

bool CBPManager::addCBP(uint i, uint j, uint k, uint l)
{
    auto key = std::make_tuple(i, j, k, l);
    if (cbp_set_.find(key) != cbp_set_.end()) {
        return false; // Already exists
    }
    
    cbp_set_.insert(key);
    cbp_.push_back(std::make_pair(std::make_pair(i, j), std::make_pair(k, l)));
    
    // Update projection arrays
    c_x_[i].push_back(j);
    c_y_[k].push_back(l);
    c_z_[i].push_back(k);
    c_z_[j].push_back(l);
    
    return true;
}

uint CBPManager::addViolatedCBPs(const std::vector<ViolationInfo>& violations)
{
    uint added = 0;
    std::set<std::tuple<uint, uint, uint, uint>> candidates_to_check;
    
    // Collect candidates based on violations
    uint total_candidates_found = 0;
    for (const auto& violation : violations) {
        auto candidates = searchCandidatesAroundViolation(violation);
        total_candidates_found += candidates.size();
        for (const auto& cand : candidates) {
            candidates_to_check.insert(cand);
        }
    }
    
    std::cout << "CBPManager: Found " << total_candidates_found << " candidates for " 
              << violations.size() << " violations, " << candidates_to_check.size() 
              << " unique candidates" << std::endl;
    
    // Check and add valid candidates
    uint rejected_by_should = 0;
    uint rejected_by_add = 0;
    for (const auto& [i, j, k, l] : candidates_to_check) {
        if (shouldAddCBP(i, j, k, l)) {
            if (addCBP(i, j, k, l)) {
                added++;
            } else {
                rejected_by_add++;
            }
        } else {
            rejected_by_should++;
        }
    }
    
    if (rejected_by_should > 0 || rejected_by_add > 0) {
        std::cout << "CBPManager: Rejected " << rejected_by_should << " by shouldAddCBP, " 
                  << rejected_by_add << " by addCBP (duplicates)" << std::endl;
    }
    
    if (added > 0) {
        std::cout << "CBPManager: Added " << added << " new CBPs based on " << violations.size() << " violations" << std::endl;
        updateProjections();
    }
    
    return added;
}

std::vector<std::tuple<uint, uint, uint, uint>> 
CBPManager::searchCandidatesAroundViolation(const ViolationInfo& violation) const
{
    std::vector<std::tuple<uint, uint, uint, uint>> candidates;
    
    if (!p_x_ || !p_y_ || !p_z_) return candidates;
    
    if (violation.type == ViolationInfo::X_VIOLATION) {
        // For X violation at (i,j), search for matching alignments
        uint i = violation.i;
        uint j = violation.j;
        
        if (i < L1_ && j < L1_ && (*p_x_)[i][j] > cutoff_ * 0.1f) { // Very relaxed cutoff
            // Search for corresponding positions in sequence 2
            for (uint k = 0; k < L2_ - 1; ++k) {
                if ((*p_z_)[i][k] <= cutoff_ * 0.1f) continue;
                
                for (uint l = k + 1; l < L2_; ++l) {
                    if ((*p_y_)[k][l] > cutoff_ * 0.1f && (*p_z_)[j][l] > cutoff_ * 0.1f) {
                        candidates.push_back(std::make_tuple(i, j, k, l));
                    }
                }
            }
        }
    }
    else if (violation.type == ViolationInfo::Y_VIOLATION) {
        // For Y violation at (k,l), search for matching alignments
        uint k = violation.k;
        uint l = violation.l;
        
        if (k < L2_ && l < L2_ && (*p_y_)[k][l] > cutoff_ * 0.1f) { // Very relaxed cutoff
            // Search for corresponding positions in sequence 1
            for (uint i = 0; i < L1_ - 1; ++i) {
                if ((*p_z_)[i][k] <= cutoff_ * 0.1f) continue;
                
                for (uint j = i + 1; j < L1_; ++j) {
                    if ((*p_x_)[i][j] > cutoff_ * 0.1f && (*p_z_)[j][l] > cutoff_ * 0.1f) {
                        candidates.push_back(std::make_tuple(i, j, k, l));
                    }
                }
            }
        }
    }
    else if (violation.type == ViolationInfo::Z_VIOLATION) {
        // For Z violation at (i,k), search for base pairs
        uint i = violation.i;
        uint k = violation.k;
        
        if (i < L1_ && k < L2_ && (*p_z_)[i][k] > cutoff_ * 0.1f) { // Very relaxed cutoff
            // Search for j paired with i in sequence 1 (forward)
            for (uint j = i + 1; j < L1_; ++j) {
                if ((*p_x_)[i][j] <= cutoff_ * 0.1f) continue;
                
                // Search for l paired with k in sequence 2
                for (uint l = k + 1; l < L2_; ++l) {
                    if ((*p_y_)[k][l] > cutoff_ * 0.1f && (*p_z_)[j][l] > cutoff_ * 0.1f) {
                        candidates.push_back(std::make_tuple(i, j, k, l));
                    }
                }
            }
            
            // Search for j paired with i in sequence 1 (backward)
            for (uint j = 0; j < i; ++j) {
                if ((*p_x_)[j][i] <= cutoff_ * 0.1f) continue;
                
                // Search for l paired with k in sequence 2
                for (uint l = k + 1; l < L2_; ++l) {
                    if ((*p_y_)[k][l] > cutoff_ * 0.1f && (*p_z_)[j][l] > cutoff_ * 0.1f) {
                        candidates.push_back(std::make_tuple(j, i, k, l));
                    }
                }
                // Search for l paired with k in sequence 2 (backward)
                for (uint l = 0; l < k; ++l) {
                    if ((*p_y_)[l][k] > cutoff_ * 0.1f && (*p_z_)[j][l] > cutoff_ * 0.1f) {
                        candidates.push_back(std::make_tuple(j, i, l, k));
                    }
                }
            }
        }
    }
    
    return candidates;
}

void CBPManager::updateProjections()
{
    // Ensure projection arrays are the correct size
    if (c_x_.size() != L1_) c_x_.resize(L1_);
    if (c_y_.size() != L2_) c_y_.resize(L2_);
    if (c_z_.size() != L1_) c_z_.resize(L1_);
    
    // Sort and remove duplicates from projection arrays
    for (uint i = 0; i < c_x_.size(); ++i) {
        std::sort(c_x_[i].begin(), c_x_[i].end());
        c_x_[i].erase(std::unique(c_x_[i].begin(), c_x_[i].end()), c_x_[i].end());
    }
    
    for (uint k = 0; k < c_y_.size(); ++k) {
        std::sort(c_y_[k].begin(), c_y_[k].end());
        c_y_[k].erase(std::unique(c_y_[k].begin(), c_y_[k].end()), c_y_[k].end());
    }
    
    for (uint i = 0; i < c_z_.size(); ++i) {
        std::sort(c_z_[i].begin(), c_z_[i].end());
        c_z_[i].erase(std::unique(c_z_[i].begin(), c_z_[i].end()), c_z_[i].end());
    }
}

bool CBPManager::cbpExists(uint i, uint j, uint k, uint l) const
{
    return cbp_set_.find(std::make_tuple(i, j, k, l)) != cbp_set_.end();
}

float CBPManager::calculateCBPScore(uint i, uint j, uint k, uint l) const
{
    if (!p_x_ || !p_y_ || !p_z_) return -1.0f;
    
    float p = (N1_ * (*p_x_)[i][j] + N2_ * (*p_y_)[k][l]) / (N1_ + N2_);
    float q = ((*p_z_)[i][k] + (*p_z_)[j][l]) / 2;
    
    return p + q; // Simplified score for sorting/prioritization
}