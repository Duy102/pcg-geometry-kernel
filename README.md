# PCG Geometry Kernel

Theorem-traceable production engineering for Prime-Curve Geometry (PCG).

## Current scope

This repository is intentionally narrow. The current implementation target is the six-arc `ABCABC` fixed-turn vertical slice from Prompt v5. It does **not** claim a complete arbitrary-trace PCG solver.

The Phase 0 audit checkpoint remains preserved in `PCG_Geometry_Kernel_Phase0_Audit.zip`.

## Phase 0–3 candidate implementation

The `phase0-3-abcabc` branch adds:

- C++20/CMake kernel skeleton;
- exact rational multiples of pi for theorem-facing turns;
- separate `InputError`, `DomainError`, `Decision`, `ProofKind`, `ArithmeticAssurance`, and `TerminationReason` semantics;
- Boost interval-based certification substrate with analytic transcendental checks;
- generic supporting-circle / finite-arc relation code and PCG `ExtraHit`;
- specialized `ABCABC` theorem solver;
- exact phase / positive-closure obstruction path;
- certified numerical cross-passage checks;
- deterministic canonical input serialization + SHA-256 binding;
- versioned certificate metadata and independent certificate re-verification;
- versioned `ABCABC` golden conformance corpus;
- GCC/Clang GitHub Actions CI.

Local clean builds were verified with GCC and Clang before publication of the candidate branch. Remote GitHub CI remains the release gate for merging.

## Scientific boundary

The stable design invariants are:

```text
REALIZABLE     => verifiable existence evidence
NOT_REALIZABLE => verifiable obstruction
otherwise      => INDETERMINATE
```

No numerical search failure is treated as proof of nonexistence. Intermediate stages must not overstate their guarantees.

## Historical Phase 0 audit

The original audit package is retained as a historical checkpoint. Its blocker statements describe the earlier environment and should not be read as the status of later implementation work.

Public visibility does not itself grant a software or manuscript license. No license has been added yet.
