# Seifert Turn-Flow Contract

## Source theorem

This module implements the relaxed feasibility layer from PCG_Seifert_Turn_Flow_v0.2.tex.

Source SHA-256:

1a60d9ea904775dfd59c345915e1cb00fb67f13dbc46a029a9a7c8846428daa5

For fixed Seifert orientations epsilon_v in {+1,-1}, piece counts P_v,N_v >= 0 with P_v+N_v >= 1, and a finite Seifert multigraph, the source defines the scaled relaxed system

    t + Bx = 2 epsilon
    -2 N_v < t_v < 2 P_v
    -1 < x_e < 1

and proves equivalence with the strict Seifert cut budgets.

## Empty-subset wording note

The v0.2 manuscript writes the strict cut inequality for every U subseteq V. Taken literally at U = empty, that becomes 0 < 0 < 0, while the Hoffman-circulation proof treats the empty cut as the vacuous 0 <= 0 condition. The implementation therefore applies the strict PCG cut test to nonempty Seifert subsets. This discrepancy is documented explicitly rather than silently changing the source theorem.

## Production decision path

The solver follows the augmented-circulation construction in the manuscript:

- non-loop Seifert edges are given one canonical arbitrary orientation;
- g_e = x_e + 1 on each non-loop edge;
- one artificial root arc r -> v is added per Seifert vertex;
- loops are discarded from net flow because their incidence contributions cancel;
- the open bounds are replaced by exact closed shrunken bounds;
- lower-bound circulation feasibility is reduced to max-flow.

Let A be the number of augmented arcs. The implementation uses the exact common margin eta = 1/(2(A+1)). The manuscript proves every satisfied strict PCG cut has integer slack at least one. Any Hoffman cut crosses at most A augmented arcs, so this shrink changes its slack by strictly less than one and preserves every strictly feasible cut.

All capacities and flows are scaled to arbitrary-precision integers. The decision path uses no floating point.

## Certificates

A FEASIBLE result carries exact scaled edge flows x_e and source-turn totals t_v with a common positive denominator D. The verifier checks directly:

    -D < x_e_scaled < D
    -2 D N_v < t_v_scaled < 2 D P_v
    t_scaled + B x_scaled = 2 D epsilon

The verifier does not rerun the max-flow solver.

An INFEASIBLE result carries a nonempty subset U and one exact violated cut:

    2 epsilon_U >= 2 P_U + d(U)

or

    2 epsilon_U <= -2 N_U - d(U)

The verifier recomputes P_U, N_U, epsilon_U and d(U) directly from the certificate-bound input.

## Validation

The CI test suite includes hand cases, the ABCADCBD-style degree-two negative leaf obstruction from v0.2, tamper tests, loop/orientation invariance, and 6000 deterministic valid-model instances with two to six vertices. For those 6000 instances the production max-flow solver is cross-checked against a separate exhaustive implementation of every nonempty strict cut inequality.

This 6000-case test is a software cross-check against an exact cut oracle. It is not claimed to reproduce the manuscript's independent LP experiment.

## Scope boundary

This is an intermediate PCG layer. FEASIBLE means only that relaxed Seifert turn totals and smoothing flows exist. It does not by itself construct individual arc turns, chord closure, prescribed Euclidean crossings, or absence of unintended intersections. It must never be promoted directly to top-level REALIZABLE.