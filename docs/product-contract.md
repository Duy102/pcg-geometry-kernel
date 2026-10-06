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
