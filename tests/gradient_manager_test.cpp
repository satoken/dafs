#include "gradient_manager.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace {

void expect(bool condition, const std::string& message)
{
    if (!condition)
        throw std::runtime_error(message);
}

void expect_same_bits(float lhs, float rhs, const std::string& message)
{
    expect(std::memcmp(&lhs, &rhs, sizeof(lhs)) == 0, message);
}

using UpdateInputs = std::tuple<
    std::vector<CBP>, VU, VU, VU, VU, VVU, VVU, VVU>;

uint apply_update(GradientManager& manager, const UpdateInputs& inputs,
                  float score = 1.0f)
{
    const auto& [cbp, x, y, z, w_cbp, c_x, c_y, c_z] = inputs;
    return manager.update_gradients(cbp, x, y, z, w_cbp, c_x, c_y, c_z,
                                    0, score);
}

UpdateInputs positive_boundary_inputs()
{
    return {{}, VU{-1u}, VU{-1u}, VU{0}, {}, VVU(1), VVU(1), VVU(1)};
}

UpdateInputs negative_z_inputs()
{
    std::vector<CBP> cbp;
    cbp.push_back({{0, 0}, {0, 0}});
    return {
        cbp,
        VU{0}, VU{0}, VU{0}, VU{0, 0},
        VVU{VU{0}}, VVU{VU{0}}, VVU{VU{0}}
    };
}

UpdateInputs positive_interior_inputs()
{
    return {{}, VU{-1u}, VU{-1u}, VU{0}, {}, VVU(1), VVU(1), VVU(1)};
}

UpdateInputs mixed_active_boundary_inputs()
{
    return {
        {},
        VU{2, -1u, -1u},
        VU{-1u, -1u, -1u},
        VU{0, 1, 2},
        {},
        VVU(3), VVU(3), VVU(3)
    };
}

void expect_nonnegative_alignment(const GradientManager& manager)
{
    expect(manager.q_z(0, 0) >= 0.0f,
           "projected Z multiplier became negative");
}

void run_projected_norm_case(bool sparse_structure, bool sparse_alignment)
{
    GradientManager manager(0.5f, 0.0f, 1.0f,
                            sparse_structure, sparse_alignment);
    manager.set_projected_norm(true);
    manager.initialize(1, 1);

    // At q_z == 0, a positive Z constraint gradient points out of the
    // feasible orthant.  It is absent from h but remains in the update.
    const UpdateInputs boundary = positive_boundary_inputs();
    const uint boundary_violations = manager.get_violations();
    apply_update(manager, boundary);
    expect(manager.get_violations() == 0,
           "positive slack at the Z boundary was counted as a violation");
    expect(manager.get_last_gradient_norm_squared() == 1.0f,
           "unexpected full norm for the boundary case");
    expect(manager.get_last_projected_gradient_norm_squared() == 0.0f,
           "outward boundary gradient was not removed from projected norm");
    expect(manager.get_last_projected_norm_dropped() == 1,
           "boundary gradient drop count is incorrect");
    expect(manager.get_last_update_size() == 0.0f &&
               std::isfinite(manager.get_last_update_size()),
           "zero projected norm did not produce a finite zero update");
    expect(manager.get_violations() == boundary_violations,
           "projected norm changed the boundary violation count");
    expect_nonnegative_alignment(manager);

    // Two identical CBP indices make the selected Z constraint negative.
    // It must remain in h even while q_z is at zero, and the projection must
    // keep the resulting multiplier feasible.
    GradientManager negative_manager(0.5f, 0.0f, 1.0f,
                                     sparse_structure, sparse_alignment);
    negative_manager.set_projected_norm(true);
    negative_manager.initialize(1, 1);
    const UpdateInputs negative = negative_z_inputs();
    apply_update(negative_manager, negative);
    expect(negative_manager.get_last_gradient_norm_squared() == 11.0f,
           "unexpected full norm for the negative Z case");
    expect(negative_manager.get_last_projected_gradient_norm_squared() ==
               11.0f,
           "negative Z gradient was removed from projected norm");
    expect(negative_manager.get_last_projected_norm_dropped() == 0,
           "negative Z case reported a dropped component");
    expect(negative_manager.get_violations() == 3,
           "projected norm changed negative/structure violation counting");
    expect(negative_manager.q_z(0, 0) > 0.0f,
           "negative Z gradient did not increase its multiplier");
    expect_nonnegative_alignment(negative_manager);
    const UpdateInputs interior = positive_interior_inputs();

    // This seed is deliberately tiny but still strictly positive.  The
    // exact q_z == 0 test must not apply an epsilon cutoff to it.
    GradientManager tiny_manager(0.5f, 0.0f, 1.0f,
                                 sparse_structure, sparse_alignment);
    tiny_manager.set_projected_norm(true);
    tiny_manager.initialize(1, 1);
    apply_update(tiny_manager, negative, 1.0e-12f);
    expect(tiny_manager.q_z(0, 0) > 0.0f &&
               tiny_manager.q_z(0, 0) < 1.0e-6f,
           "strictly positive q_z was not preserved at tiny scale");
    apply_update(tiny_manager, interior);
    expect(tiny_manager.get_last_gradient_norm_squared() == 1.0f &&
               tiny_manager.get_last_projected_gradient_norm_squared() ==
                   1.0f &&
               tiny_manager.get_last_projected_norm_dropped() == 0,
           "tiny positive q_z was incorrectly treated as a boundary");
    expect_nonnegative_alignment(tiny_manager);

    // The same positive gradient is now at q_z > 0, so it remains active in
    // h and is not dropped merely because it is positive.
    apply_update(negative_manager, interior);
    expect(negative_manager.get_last_gradient_norm_squared() == 1.0f,
           "unexpected full norm for the positive interior case");
    expect(negative_manager.get_last_projected_gradient_norm_squared() ==
               1.0f,
           "positive gradient at q_z > 0 was incorrectly removed");
    expect(negative_manager.get_last_projected_norm_dropped() == 0,
           "positive interior case reported a dropped component");
    expect(negative_manager.get_violations() == 0,
           "positive interior case changed violation counting");
    expect_nonnegative_alignment(negative_manager);
}

void run_mixed_active_boundary_case(bool sparse_structure,
                                    bool sparse_alignment)
{
    const UpdateInputs inputs = mixed_active_boundary_inputs();

    GradientManager projected(0.5f, 0.0f, 1.0f,
                              sparse_structure, sparse_alignment);
    projected.set_projected_norm(true);
    projected.initialize(3, 3);
    const uint projected_violations = apply_update(projected, inputs);
    expect(projected_violations == 1,
           "mixed fixture changed its projected violation count");
    expect(projected.get_last_gradient_norm_squared() == 4.0f &&
               projected.get_last_projected_gradient_norm_squared() == 1.0f,
           "mixed fixture did not isolate the active X norm");
    expect(projected.get_last_projected_norm_dropped() == 3,
           "mixed fixture dropped the wrong number of boundary Z entries");
    expect(projected.get_last_update_size() == 0.5f,
           "projected mixed fixture has an unexpected Polyak step");
    expect(projected.q_x(0, 2) == 0.5f,
           "projected mixed fixture did not use the active norm step");
    for (uint k = 0; k != 3; ++k)
        expect(projected.q_z(0, k) == 0.0f,
               "outward boundary Z multiplier changed in mixed fixture");
    expect_nonnegative_alignment(projected);

    // Fixed values from the legacy path make this a regression reference,
    // rather than only a comparison between two new flag settings.
    GradientManager legacy(0.5f, 0.0f, 1.0f,
                           sparse_structure, sparse_alignment);
    legacy.initialize(3, 3);
    expect(!legacy.uses_projected_norm(),
           "legacy reference unexpectedly enabled projected norm");
    expect(apply_update(legacy, inputs) == 1,
           "legacy mixed fixture changed its violation count");
    expect(legacy.get_last_gradient_norm_squared() == 4.0f &&
               legacy.get_last_projected_gradient_norm_squared() == 4.0f &&
               legacy.get_last_projected_norm_dropped() == 0,
           "legacy mixed fixture norm reference changed");
    expect(legacy.get_last_update_size() == 0.125f &&
               legacy.q_x(0, 2) == 0.125f,
           "legacy mixed fixture Polyak trajectory changed");
    for (uint k = 0; k != 3; ++k)
        expect(legacy.q_z(0, k) == 0.0f,
               "legacy boundary Z multiplier changed in mixed fixture");
}

void run_legacy_disabled_case(bool sparse_structure, bool sparse_alignment)
{
    GradientManager implicit_false(0.5f, 0.0f, 1.0f,
                                   sparse_structure, sparse_alignment);
    GradientManager explicit_false(0.5f, 0.0f, 1.0f,
                                   sparse_structure, sparse_alignment);
    explicit_false.set_projected_norm(false);
    implicit_false.initialize(1, 1);
    explicit_false.initialize(1, 1);

    const std::vector<UpdateInputs> sequence = {
        positive_boundary_inputs(), negative_z_inputs(),
        positive_interior_inputs()
    };
    for (const UpdateInputs& inputs : sequence) {
        apply_update(implicit_false, inputs);
        apply_update(explicit_false, inputs);

        expect(!implicit_false.uses_projected_norm() &&
                   !explicit_false.uses_projected_norm(),
               "legacy managers unexpectedly enabled projected norm");
        expect(implicit_false.get_violations() ==
                   explicit_false.get_violations(),
               "legacy violation count changed between equivalent modes");
        expect_same_bits(
            implicit_false.get_last_update_size(),
            explicit_false.get_last_update_size(),
            "legacy Polyak update is not bitwise identical");
        expect_same_bits(
            implicit_false.get_last_gradient_norm_squared(),
            explicit_false.get_last_gradient_norm_squared(),
            "legacy full norm is not bitwise identical");
        for (uint i = 0; i != 1; ++i) {
            for (uint j = 0; j != 1; ++j) {
                expect_same_bits(implicit_false.q_x(i, j),
                                 explicit_false.q_x(i, j),
                                 "legacy q_x changed with the flag disabled");
                expect_same_bits(implicit_false.q_y(i, j),
                                 explicit_false.q_y(i, j),
                                 "legacy q_y changed with the flag disabled");
                expect_same_bits(implicit_false.q_z(i, j),
                                 explicit_false.q_z(i, j),
                                 "legacy q_z changed with the flag disabled");
            }
        }
    }
}

void run_sparse_raw_visitor_case()
{
    SparseFloatMatrix matrix;
    matrix.assign(2, 1000000);
    matrix.set(0, 7, 3.0f);
    matrix.set(0, 999999, 2.0f);
    matrix.set(0, 500003, 4.0f);
    matrix.set(1, 999998, -1.0f);

    size_t count = 0;
    double sum = 0.0;
    bool saw_concentrated_row = false;
    matrix.for_each_nonzero([&](uint row, uint column, float value) {
        expect(row < 2 && column < 1000000,
               "raw sparse visitor returned an invalid coordinate");
        ++count;
        sum += static_cast<double>(value);
        if (row == 0 && (column == 7 || column == 500003 ||
                         column == 999999))
            saw_concentrated_row = true;
    });
    expect(count == 4, "raw sparse visitor missed a concentrated entry");
    expect(sum == 8.0, "raw sparse visitor returned an incorrect sum");
    expect(saw_concentrated_row,
           "raw sparse visitor did not cover the concentrated row");

    matrix.set(0, 7, 0.0f);
    count = 0;
    matrix.for_each_nonzero([&](uint, uint, float) { ++count; });
    expect(count == 3, "raw sparse visitor retained a removed entry");
}

} // namespace

int main()
{
    try {
        run_sparse_raw_visitor_case();
#if defined(USE_ADAGRAD) || defined(USE_ADAM)
        GradientManager manager(0.5f);
        bool rejected = false;
        try {
            manager.set_projected_norm(true);
        } catch (const std::invalid_argument&) {
            rejected = true;
        }
        expect(rejected,
               "adaptive gradient mode accepted unsupported projected norm");
#else
        for (const bool sparse_structure : {false, true}) {
            for (const bool sparse_alignment : {false, true}) {
                run_projected_norm_case(sparse_structure, sparse_alignment);
                run_mixed_active_boundary_case(sparse_structure,
                                               sparse_alignment);
                run_legacy_disabled_case(sparse_structure, sparse_alignment);
            }
        }
#endif
    } catch (const std::exception& error) {
        std::cerr << "gradient_manager_test: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
