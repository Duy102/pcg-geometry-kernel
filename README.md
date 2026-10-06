# PCG Geometry Kernel

Theorem-traceable production engineering for Prime-Curve Geometry (PCG).

## Current status

The repository currently contains a **hardened certified ABCABC vertical slice**, an exact **Seifert turn-flow** layer, a **general certified Network Closure** module, a **Phase 6A projectively-rigid post-closure geometry** verifier, and a **Phase 6B1 higher-dimensional positive-kernel witness** path, and a **Phase 6B2A solver-neutral semialgebraic IR compiler**.

It is intentionally narrow: this is **not** a claim that arbitrary-trace PCG or the full PCG framework is production-complete.

The historical Phase 0 audit checkpoint remains preserved in `PCG_Geometry_Kernel_Phase0_Audit.zip`.

## Implemented ABCABC vertical slice

The current `main` branch includes:

- C++20/CMake kernel;
- arbitrary-precision exact rational multiples of pi for theorem-facing turns;
- explicit `InputError`, `DomainError`, `Decision`, `ProofKind`, `ArithmeticAssurance`, and `TerminationReason` semantics;
- Boost interval-based certification substrate;
- generic supporting-circle / finite-arc relation code and PCG `ExtraHit`;
- specialized fixed-turn `ABCABC` theorem solver;
- exact paired phase and positive-chord-closure obstruction paths;
- certified numerical cross-passage checks;
- explicit theorem genericity/domain gate for required support-circle distinctness, prescribed transversality, and antipodal branch separation;
- deterministic canonical input serialization and SHA-256 certificate binding;
- theorem/source provenance metadata;
- certificate verifier in a separate translation unit from the solver's private evaluator;
- tamper checks for theorem identity, source digest, input digest, proof reason, proof kind, and arithmetic assurance;
- versioned fixed-turn golden conformance corpus;
- exhaustive `ABCABC` sign-rotation conformance suite.

## Exhaustive classification coverage

The sign-rotation theorem suite covers all:

```text
64 sign patterns x 3 rotations {-2, 0, 2} = 192 classes
```

The expected realizable counts are:

- `r = -2`: 15;
- `r = 0`: 6;
- `r = 2`: 15;
- total: 36.

The theorem's rational orbit witnesses are expanded by cyclic shift and global sign reflection so that all 36 realizable sign-rotation classes are exercised through the fixed-turn solver and certificate verifier.

## Scientific decision contract

```text
REALIZABLE     => theorem-domain assumptions certified + verifiable existence evidence
NOT_REALIZABLE => verifiable theorem obstruction
UNSUPPORTED    => certified violation of a required theorem-domain assumption
INDETERMINATE  => certification cannot separate a required numerical/geometric boundary
```

No numerical search failure is treated as proof of nonexistence.

## Phase 4 Seifert turn-flow layer

The production-candidate stack now also contains an exact relaxed Seifert turn-flow module based on PCG Seifert Turn-Flow Theory v0.2.

It uses arbitrary-precision integer max-flow after an exact lower-bound circulation reduction. FEASIBLE certificates contain exact interior flow/turn witnesses; INFEASIBLE certificates contain exact violated nonempty Seifert cuts. The verifier checks these witnesses algebraically without rerunning max-flow.

The module is intentionally intermediate: Seifert-flow feasibility does not imply full geometric realizability.

Its validation includes a deterministic 6000-case cross-check against a separate exhaustive exact-cut oracle.

## Phase 5A–5B Network Closure

For rational-pi turns, the kernel builds the Network Closure matrix with exact cyclotomic algebra rather than ordinary floating-point sine/cosine values.

Phase 5A introduced exact matrix rank/nullity, exact full-rank obstruction, and certified projectively-rigid (`dim ker N = 1`) sign decisions.

Phase 5B extends this to higher-dimensional kernels. A closed instance carries an explicit algebraic positive-kernel witness `c > 0` with `N c = 0`; a non-closed instance carries an explicit Stiemke dual witness `y` with `N^T y >= 0` and `N^T y != 0`. Certificate schema v1.1 verifies these witnesses directly.

Validation includes the production ABCABC witness plus regular and open-half-plane cycle families exercising nullity 2 through 6.

The deliberate algebraic software boundary remains cyclotomic order `<= 256`; unresolved real-sign separation returns `INDETERMINATE` rather than a guessed decision.

## Phase 6A projectively-rigid post-closure geometry

Phase 6A implements the first production specialization of the PCG General Fixed-Turn Intersection Feasibility theorem.

For a certified closed network with `dim ker N = 1`, the kernel reconstructs the unique normalized metric skeleton in exact cyclotomic coordinates. It then verifies:

- prescribed-passage multiplicity and tangent transversality;
- pairwise distinct quotient vertices;
- every unordered finite circular-arc pair;
- only the quotient-vertex intersections prescribed by the trace are allowed.

Remote pairs use a general certified finite-arc `Hit` predicate. One-endpoint pairs use `ExtraHit`, with an exact tangent-line fast path for distinct support circles. Finite-arc membership follows the theorem's minor / semicircle / major cross-product formulas.

Tangency, coincident-support, or interval boundaries that cannot be certified return `INDETERMINATE`. Higher-dimensional post-closure length selection is deliberately left for a later phase because the trace-faithful subset of the positive closure polytope can be nonconvex.

## Phase 6B1 higher-dimensional positive-kernel witnesses

For `dim ker N > 1`, Phase 5B already provides one exact algebraic positive-kernel witness. Phase 6B1 can now reconstruct that specific normalized metric skeleton and run the same complete finite-arc pair verifier used by Phase 6A.

The logic is deliberately asymmetric:

- if the supplied positive-kernel metric passes every trace-faithfulness check, it is a verifiable **REALIZABLE** existence certificate;
- if that metric has a vertex collision or unintended intersection, the result is only **INDETERMINATE** for the global problem, because another positive kernel vector may still work;
- an upstream Stiemke / full-rank Network Closure obstruction remains a global **NOT_REALIZABLE** certificate.

This implements a sound higher-dimensional existence path without pretending that one sampled point solves the theorem's generally nonconvex length-selection problem. Full logical completeness for higher-dimensional fixed-turn traces still requires a finite semialgebraic decision backend, e.g. real quantifier elimination / CAD, or a mathematically equivalent specialized solver.

## Phase 6B2A semialgebraic backend contract

The kernel now compiles every rational-pi fixed-turn trace into a deterministic, versioned solver-neutral IR matching the General Fixed-Turn Intersection Feasibility theorem's quantified sentence.

The IR records exact tangent/chord phases, minor/semicircle/major arc branches, all quotient-vertex distinctness constraints, every unordered bad-pair incidence formula with its allowed common vertices, prescribed fixed-tangent transversality checks, and quantifier counts. It is bound to the canonical input and theorem source by SHA-256 and can be independently recompiled and structurally verified.

This phase freezes the theorem-to-solver boundary. It does **not** yet trust or embed a CAD / real-quantifier-elimination engine, and therefore does not add new global decisions beyond the existing certified paths.

## CI / Trusted Computing Base gates

The current CI baseline uses:

- Ubuntu 24.04;
- Boost 1.83 package `libboost1.83-dev=1.83.0-2.1ubuntu3.2`;
- GCC 13 and Clang 18;
- a pinned immutable `actions/checkout` commit.

Every push / pull request runs:

- `linux-gcc-Debug`;
- `linux-gcc-Release`;
- `linux-clang-Debug`;
- `linux-clang-Release`;
- `linux-gcc-asan-ubsan` (AddressSanitizer + UndefinedBehaviorSanitizer).

The hardened ABCABC, Seifert, and general Network Closure merges passed all five gates. Phase 6A is held to the same five-gate requirement before merge.

This improves reproducibility but is not a fully hermetic build: the GitHub Ubuntu image and package repository can still evolve.

## Deliberate remaining scope boundaries

The current repository does not yet provide:

- an evidence-verifiable exact backend for the Phase 6B2A semialgebraic IR (real quantifier elimination / CAD or an equivalent exact solver);
- production `ABCADCBD` specialization;
- global BVH / spatial-index acceleration for large remote-pair workloads;
- variable-turn / sign-rotation general decision machinery;
- language bindings;
- a stable public v1 API.

These belong to later PCG phases rather than the present ABCABC vertical slice.

## Licensing

Public visibility does not itself grant a software or manuscript license. No license has been added yet; licensing remains an explicit project decision.
