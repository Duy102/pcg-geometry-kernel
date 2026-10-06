# Phase 0–3 Status

Implemented under Prompt v5:
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

Publication note: the reviewed Phase 0–3 source archive was reconstructed on GitHub and verified against SHA-256 `9c02c915fd689b60640531710c0c690ff40dfd95dbb46a4333d532973a9bad61` before expansion.

Current deliberate limitation: theorem-facing input uses exact rational multiples of pi. This avoids false exact-modulo claims for decimal-radian approximations. Broader real-input encodings are future work.

Remote GCC/Clang CI is the merge gate for this branch.
