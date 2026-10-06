# Phase 6A — Projectively-Rigid Post-Closure Geometry

## Theorem basis

Source: `PCG_General_Fixed_Turn_Intersection_Feasibility_v0.1.tex`

SHA-256:

`625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208`

The source theorem states that after positive Network Closure, trace-faithfulness
is obtained by requiring distinct quotient vertices and excluding every
unintended common point of finite circular-arc pairs. Finite-arc membership is
semialgebraic:

- minor arc: `A >= 0 && B >= 0`;
- semicircle: `A >= 0`;
- major arc: `A >= 0 || B >= 0`.

For rational multiples of pi the fixed trigonometric data are algebraic.

## Phase 6A specialization

This implementation deliberately starts with the projectively rigid case:

`dim ker N = 1`.

After scale normalization, the positive Network Closure certificate determines
one exact algebraic metric skeleton. The kernel exports that skeleton in the
same cyclotomic power basis used by the Network Closure certificate. Geometry
is then evaluated with directed-rounding intervals; ordinary floating-point
coordinates are never promoted to proof evidence.

The verifier checks:

1. the upstream Network Closure certificate;
2. prescribed-passage multiplicity/transversality;
3. pairwise distinct quotient vertices;
4. every unordered finite-arc pair.

For pairs with one prescribed common endpoint, the existing PCG `ExtraHit`
predicate is used. For pairs with two prescribed endpoints, distinct support
circles imply that the two allowed endpoints exhaust their circle
intersections. For remote pairs, Phase 6A adds a general certified finite-arc
`Hit` predicate based on circle intersections plus the theorem's
minor/semicircle/major cross-product membership formulas.

Tangency, coincident-support, or interval-separation boundaries that are not
certified are returned as `INDETERMINATE`; they are never guessed.

## Boundary

This is not yet the full higher-dimensional theorem implementation. When
`dim ker N > 1`, trace-faithful length selection can be nonconvex; Phase 6A
therefore returns `HigherDimensionalClosure / INDETERMINATE`. A future Phase
6B may add semialgebraic/low-dimensional search or a more specialized
certificate scheme.
