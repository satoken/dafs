# Opt-in tight DD alignment bound

`--dd-tight-alignment` enables an exact sparse oracle for the positive part of
the monotone alignment subproblem. It is off by default. With the option off,
the existing row/column maxima, certified-track row relaxation, repair fields,
and metric key sets are unchanged.

For a sparse support edge `(i,k)`, the oracle uses

```text
w(i,k) = p_z(i,k) - align-th + q_z(i,k)
```

and maximizes the sum of positive weights over matchings whose rows and
columns are both strictly increasing. Edges are processed in row groups. A
Fenwick tree stores the best predecessor ending in a column strictly smaller
than the current column, giving `O(m log L)` time and `O(m + L)` storage for
`m` supported edges and profile dimensions bounded by `L`. Delaying all updates
from one row prevents two edges from the same row being chained together.
Coordinates are canonicalized and ties are deterministic; the returned
mapping is also an exact subgradient of this maximum-of-affine function.

The opt-in support is the union of nonzero probability entries and nonzero
`q_z` entries from the active and certified tracks. This includes q-only
edges. If an edge is absent from both sources, its score is `-align-th`; the
option therefore requires a finite, non-negative `--align-th`, so omitted
edges cannot have positive weight. A negative or non-finite threshold is
rejected during option parsing rather than producing an uncertified bound.
Scores are accumulated in `double` and rounded outward to `float` before
being used as an upper bound.

The same switch tightens the main alignment upper bound and replaces the
certified track’s row-capacity alignment oracle with the exact monotone
oracle. The certified track’s structure components intentionally remain the
existing left-endpoint relaxation. It is therefore still possible to have a
nonzero duality gap, and beam decoding is not made globally exact. Exact CBP
pricing remains in place: q-only alignment edges participate in the matching
oracle, but cannot form a valid CBP under the current candidate validator
without the required positive probability endpoints; `q_z` also enters
reduced cost with the sign that cannot create a missing positive column by
itself.

The monotone certified mapping is tested as an additional feasible repair. Its
score is reported as `certified_alignment_repair`; the existing
`feasible_repair` field remains exactly `max(intersection_repair,
consensus_repair)`. Its repair time is included in DD time and is reported in
the opt-in metrics. The feasible primal lower bound may use the additional
candidate.

This is an acceleration/certificate tightening, not a convergence theorem.
The current Polyak update uses a feasible primal lower bound, but integrality
and relaxation gaps can prevent dual-optimum convergence. Fewer DD iterations
are an empirical possibility, not a guarantee; total runtime, including the
extra repair, remains the measured quantity. Because the option can change
the bound, lower bound, dual trajectory, and final prediction, it is isolated
behind the opt-in switch. For the underlying theory, see Stephen Boyd and
Jaehyun Park, [*Subgradient Methods*](https://stanford.edu/class/ee364b/lectures/subgrad_method_notes.pdf),
especially §§4 and 6 on Polyak step lengths and projected subgradient methods
for dual problems. In this implementation the lower-bound target is a feasible
primal value, not necessarily the dual optimum; the existing constant Polyak
scaling therefore does not by itself provide a convergence guarantee.
