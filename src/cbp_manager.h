/*
 * CBP Manager for Column Generation in DAFS
 * Manages dynamic addition of consensus base pairs based on constraint violations
 */

#ifndef CBP_MANAGER_H
#define CBP_MANAGER_H

#include <vector>
#include <set>
#include <tuple>
#include <unordered_map>
#include <cstddef>
#include "typedefs.h"
#include "gradient_manager.h"

class CBPManager {
public:
    using CBP = std::pair<std::pair<uint, uint>, std::pair<uint, uint>>;
    using VF = std::vector<float>;
    using VVF = std::vector<VF>;
    using VU = std::vector<uint>;
    using VVU = std::vector<VU>;
    
    // ViolationInfo is defined in gradient_manager.h

private:
    std::vector<CBP> cbp_;
    std::set<std::tuple<uint, uint, uint, uint>> cbp_set_;
    VVU c_x_, c_y_, c_z_;
    
    // Problem dimensions
    uint L1_, L2_;
    
    // Thresholds and parameters
    float cutoff_;
    float w_;
    float th_a_;
    std::vector<float> th_s_;
    
    // Reference to probability matrices
    const VVF* p_x_;
    const VVF* p_y_;
    const VVF* p_z_;
    uint N1_, N2_;

public:
    CBPManager(uint L1, uint L2, float cutoff = 0.01f);
    
    // Initialize with probability matrices and parameters
    void initialize(const VVF& p_x, const VVF& p_y, const VVF& p_z,
                   uint N1, uint N2, float w, float th_a, const std::vector<float>& th_s);
    
    // Generate initial CBP set with strict thresholds
    void generateInitialCBPs();
    
    // Add new CBPs based on violations
    uint addViolatedCBPs(const std::vector<ViolationInfo>& violations);
    
    // Check if a CBP candidate should be added
    bool shouldAddCBP(uint i, uint j, uint k, uint l) const;
    
    // Add a single CBP if not already present
    bool addCBP(uint i, uint j, uint k, uint l);
    
    // Getters
    const std::vector<CBP>& getCBPs() const { return cbp_; }
    const VVU& get_c_x() const { return c_x_; }
    const VVU& get_c_y() const { return c_y_; }
    const VVU& get_c_z() const { return c_z_; }
    size_t size() const { return cbp_.size(); }
    
    // Search for violated CBP candidates around a given position
    std::vector<std::tuple<uint, uint, uint, uint>> 
    searchCandidatesAroundViolation(const ViolationInfo& violation) const;
    
    // Update projection arrays after adding CBPs
    void updateProjections();
    
private:
    // Helper to check if CBP exists
    bool cbpExists(uint i, uint j, uint k, uint l) const;
    
    // Calculate CBP score
    float calculateCBPScore(uint i, uint j, uint k, uint l) const;
};

#endif // CBP_MANAGER_H