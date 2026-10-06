# Phase 6B1 — Higher-Dimensional Positive-Kernel Witness Path

## Theorem basis

Source: `PCG_General_Fixed_Turn_Intersection_Feasibility_v0.1.tex`

SHA-256:

`625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208`

The general theorem identifies trace-faithful fixed-turn realizations with
points `(p,c)` in a semialgebraic feasible set. For rational multiples of
pi, all required fixed trigonometric coefficients are algebraic and the full
decision problem is exactly decidable by real quantifier elimination.

For `dim ker N > 1`, the normalized positive closure polytope has positive
dimension and the forbidden Hit regions need not be convex. Therefore a failed
metric witness is not a global non-realizability proof.

## Production step in 6B1

Phase 5B Network Closure already emits one exact algebraic vector

`c > 0, N c = 0`.

Phase 6B1 now:

1. verifies that upstream certificate;
2. normalizes its exact positive-kernel witness by `sum(c)=1`;
3. reconstructs quotient vertices in exact cyclotomic coordinates;
4. certifies prescribed crossing transversality and distinct vertices;
5. checks every unordered finite-arc pair with `Hit` / `ExtraHit`;
6. emits a theorem/source-bound post-closure witness certificate.

If all checks pass, this one point belongs to the theorem's trace-faithful set,
so `REALIZABLE` is certified.

If the chosen point has a collision or unintended intersection, Phase 6B1
returns `INDETERMINATE`, not `NOT_REALIZABLE`. The positive closure
polytope may contain another trace-faithful point.

## Remaining completeness gap

The theorem's full higher-dimensional decision sentence is schematically

`exists (p,c) [Network && Positive && Normalize && DistinctVertices && all pair exclusions]`.

Closing this gap requires an exact semialgebraic decision backend such as real
quantifier elimination / cylindrical algebraic decomposition, or a
mathematically equivalent specialized solver. Phase 6B1 intentionally does
not claim that one supplied witness replaces that quantified search.
