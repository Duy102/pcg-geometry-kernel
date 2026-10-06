# PCG Geometry Kernel

Theorem-traceable production engineering for Prime-Curve Geometry (PCG).

## Current status

The repository currently contains a **hardened certified ABCABC production-candidate vertical slice** merged into `main`.

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

The hardened ABCABC merge and the post-merge `main` run passed all five gates.

This improves reproducibility but is not a fully hermetic build: the GitHub Ubuntu image and package repository can still evolve.

## Deliberate remaining scope boundaries

The current repository does not yet provide:

- a general arbitrary-Gauss-word solver;
- the full general Seifert / Network Closure stack;
- production `ABCADCBD`;
- general remote-Hit / global BVH or spatial-index infrastructure;
- language bindings;
- a stable public v1 API.

These belong to later PCG phases rather than the present ABCABC vertical slice.

## Licensing

Public visibility does not itself grant a software or manuscript license. No license has been added yet; licensing remains an explicit project decision.
