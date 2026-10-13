# Opt-in DD projected-active norm

`--dd-projected-norm` is a low-cost, opt-in DD step-size pilot. It changes
only the Polyak denominator; it does not change either decoder, exact CBP
pricing, the stopping tolerance, violation counting, clipping, or the
projected multiplier update. With the option absent, the legacy norm and
metric key set are unchanged.

For the materialized subgradient `g`, define `h` coordinatewise by

```text
h_z(i,k) = 0  if q_z(i,k) == 0 and g_z(i,k) > 0
           g_z(i,k) otherwise
h_x       = g_x
h_y       = g_y
```

The equality test is exact: small positive multipliers are never treated as
zero. The Polyak denominator is `||h||^2` when the option is enabled. The
update still uses the existing `g`, clipping, and `q_z = max(0, q_z -
update)`, so no outward-at-boundary coordinate can make `q_z` negative. At a
boundary coordinate, projecting a step with `g_z > 0` or with the dropped
`h_z = 0` gives the same zero coordinate for every nonnegative step length.

## Conditional theory

Let `f` be the convex minimization form of the dual objective, let `q_z >= 0`,
and let `v` be any feasible point. If `g` satisfies the subgradient inequality,
then every removed coordinate contributes

```text
g_z(i,k) * (v_z(i,k) - q_z(i,k)) >= 0
```

because it is removed only at `q_z(i,k) = 0` and every feasible `v_z(i,k)` is
nonnegative. Therefore `h` also satisfies the relative-domain inequality
`f(v) >= f(q) + h . (v-q)`. This is the projected-subgradient rationale for
using `||h||^2` in the denominator; it is not an assertion that the beam
direction is an exact subgradient of an uncomputed oracle.

The standard projected-subgradient/Polyak argument applies conditionally when
the objective, subgradient, target optimum, and un-clipped update assumptions
required by that argument hold. This implementation uses a feasible primal
lower bound rather than a known dual optimum, keeps update clipping, and may
use beam directions. Consequently this option preserves bound and feasibility
correctness but does not provide a general convergence theorem or guarantee
fewer iterations. See Boyd and Park, *Subgradient Methods*, §§4 and 6:
[Stanford EE364b notes](https://stanford.edu/class/ee364b/lectures/subgrad_method_notes.pdf).

## Diagnostics and validation

When enabled, each `dd_iteration` event reports the full `g^2`, projected
`h^2`, and exact dropped-entry count for the beam track. The same three fields
are reported for the certified track when that track exists:

```text
beam_gradient_norm_squared
beam_projected_gradient_norm_squared
beam_projected_norm_dropped
certified_gradient_norm_squared
certified_projected_gradient_norm_squared
certified_projected_norm_dropped
```

`tests/gradient_manager_test.cpp` covers dense and sparse structure/alignment
storage, exact-zero boundary drops, mixed active/boundary gradients, positive
and negative Z gradients, tiny-but-positive multipliers, finite zero-norm
steps, nonnegative Z multipliers, unchanged violations, and a frozen
one-iteration legacy trajectory. It also checks that implicit and explicit
disabled settings are bitwise equivalent. Adaptive compile-time update modes
reject the option because their update does not use the Polyak denominator.

This is intentionally the first pilot. The prior `--dd-tight-alignment`
oracle remains separate and is not folded into this candidate; no periodic
tight-bound evaluation is enabled by this change.

## Opt-in certified column-capacity bound

`--dd-column-bound` is a separate, norm-off pilot for `LinearAlign`. It
requires a finite, non-negative `--align-th` and is rejected for nonlinear or
auxiliary alignment decoders. The option does not change the row-capacity
certified track, its selected alignment, Polyak value, gradient, pricing, or
lower bound. It only adds a second valid prefix upper bound to `best_ub` and
therefore can affect certified-gap stopping.

At the current certified `q`, the extra bound is accumulated as

```text
structure-left-X + structure-left-Y
+ column-capacity-alignment(p_z union certified nonzero q_z)
+ the same exactly priced positive CBP terms.
```

All terms use double accumulation and outward float rounding. A monotone
alignment is a matching, so its value is at most the sum of the largest
nonnegative edge in each column. The support includes q-only edges; absent
edges have score `-align-th <= 0`, which is why the threshold precondition is
required. Duplicate support entries are harmless for a maximum and are
appended without a per-iteration hash set. The existing row-based bound and,
when selected, the tight monotone bound remain independent candidates; the
minimum is used only for the prefix certificate. The option is disabled by
default and is not combined with `--dd-projected-norm`.
