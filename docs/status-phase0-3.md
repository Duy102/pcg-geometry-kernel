# Phase 0–3 Status

Implemented locally under Prompt v5:
- C++20/CMake production skeleton;
- exact rational multiples of pi for theorem-facing turns;
- `InputError`/`DomainError`/Decision/ProofKind/ArithmeticAssurance/TerminationReason separation;
- Boost interval certification substrate and analytic sin/cos tests;
- finite shared-start circle relation and PCG `ExtraHit` layer;
- specialized `ABCABC` theorem solver;
- exact phase/positive-closure obstruction path;
- certified numerical cross-passage path;
- versioned certificate metadata and theorem/source digest;
- canonical input serialization + SHA-256 binding;
- verifier that recomputes the theorem result rather than trusting a stored status;
- versioned golden conformance corpus.

Current deliberate limitation: theorem-facing input uses exact rational multiples of pi. This avoids false exact-modulo claims for decimal-radian approximations. Broader real-input encodings are future work.
