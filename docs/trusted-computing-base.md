# Trusted Computing Base

For Phase 0–3 certified results, the verifier relies on:
- C++20 control flow plus Boost.Multiprecision `cpp_int` for arbitrary-precision integer bookkeeping;
- Boost.Rational over `cpp_int` for theorem-facing exact rational multiples of pi;
- Boost.Multiprecision `cpp_dec_float<100>` only as an intermediate when widening exact rationals into binary64 certification intervals;
- Boost.Numeric.Interval;
- platform floating-point rounding control (`fenv`) selected by Boost;
- platform libm behavior used through Boost `rounded_transc_std`;
- explicit pi bracketing constants;
- certificate parser/serializer once external serialization is enabled;
- core finite-arc relation implementation used by the verifier.

Differential tests validate logic above this shared numerical substrate; they do not independently prove the TCB itself. Analytic known-value tests are therefore required for the TCB.

## Certificate verifier separation

The production ABCABC certificate verifier is compiled from a separate translation unit and does not call the solver's private `evaluate()` path. It independently recomputes phase congruences, positive closure, and the six cross-passage decisions, and it rejects tampering of the stored proof reason, proof kind, or arithmetic assurance. The verifier deliberately still shares the documented exact-arithmetic, interval, transcendental, and shared-start geometry TCB; this is implementation-path separation, not a claim of a fully independent numerical backend.

## CI reproducibility baseline

The merge-gate workflow is pinned to `ubuntu-24.04`, Boost headers package `libboost1.83-dev=1.83.0-2.1ubuntu3.2`, compiler majors GCC 13 / Clang 18, and an immutable `actions/checkout` commit SHA. Every matrix run prints a toolchain manifest containing compiler, CMake, Boost package, OpenSSL, and glibc versions. Both Debug and Release are tested on GCC and Clang. A separate GCC Debug job runs AddressSanitizer plus UndefinedBehaviorSanitizer (ASan/UBSan).

This materially reduces CI drift but is not a claim of bit-for-bit hermetic builds: GitHub's `ubuntu-24.04` runner image and Ubuntu package repository can still evolve. A container-by-digest or archived package snapshot would be the next step if full hermetic reproducibility becomes a release requirement.
