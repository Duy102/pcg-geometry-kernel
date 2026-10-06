# Product Contract

## v1 direction
PCG Geometry Kernel provides general fixed-turn PCG primitives and certified end-to-end solvers for selected theorem-supported finite circular-arc trace families.

## Immediate supported solver target
`ABCABC` under the assumptions of `PCG Complete Six-Arc ABCABC Realizability v0.3`:
- six source-ordered circular arcs with Gauss word `ABCABC`;
- `0 < |tau_i| < 2*pi`;
- C1 source passage;
- prescribed double points transverse;
- required distinct supporting circles are distinct;
- no antipodal/tangent degeneracy except prescribed source tangencies.

The initial certified input encoding represents each turn as an exact rational multiple of pi. This is a software-domain restriction, not a restriction of the theorem.

## Global invariants
- `REALIZABLE` implies verifiable existence evidence.
- `NOT_REALIZABLE` implies a verifiable obstruction.
- otherwise a valid supported instance is `INDETERMINATE`.

Intermediate modules never emit top-level realizability decisions.

## ABCABC theorem-domain gate

After phase and positive-chord feasibility, the specialized solver explicitly certifies the v0.3 genericity assumptions before any top-level REALIZABLE claim. The gate checks required support-circle distinctness across the fifteen arc pairs, cross-passage transversality at prescribed double points, and separation from the antipodal principal-argument branch cut. Certified violations return UNSUPPORTED; unresolved interval separation returns INDETERMINATE.

## Certificate integrity

Verification does not trust certificate claim fields. The verifier recomputes the supported ABCABC theorem path independently from the solver's private evaluator and requires the stored reason, proof kind, arithmetic assurance, theorem identity, source digest, and canonical input digest to agree with recomputation. Shared low-level geometry and certified numerical primitives remain part of the common TCB.

## Relaxed Seifert turn-flow layer

The Phase 4 Seifert module is an intermediate certified feasibility layer, not a top-level geometric realizability solver.

For a finite Seifert multigraph with fixed circle orientations and valid positive/negative source-piece capacities, it returns:

- an exact interior turn-flow witness when the relaxed system is feasible;
- an exact nonempty Seifert cut obstruction when the relaxed system is infeasible;
- INDETERMINATE only for backend or certificate-integrity failure.

The module must never map relaxed flow feasibility directly to top-level REALIZABLE. Individual arc-turn allocation, metric closure, prescribed crossing geometry, and unintended-intersection exclusion remain later layers.

## Network Closure Phase 5A

The Network Closure module is an intermediate metric-closure layer. It consumes a cyclic quotient trace and exact rational-pi turns and tests the source theorem condition `ker(N) intersect R^m_{>0}` for the currently supported exact algebraic cases.

Phase 5A is complete for:

- exact full-column-rank obstruction;
- projectively rigid networks (`dim ker N = 1`) when algebraic signs are certifiably separated.

It returns `INDETERMINATE` for higher-dimensional kernels, cyclotomic order above the software cap, or unresolved sign separation. In particular, metric `CLOSED` must never be promoted directly to top-level `REALIZABLE`; finite-arc collision/intersection checks remain separate.
