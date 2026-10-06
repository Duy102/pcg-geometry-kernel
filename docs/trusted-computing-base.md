# Trusted Computing Base

For Phase 0–3 certified results, the verifier relies on:
- C++20 integer/rational arithmetic used for exact pi-multiple bookkeeping;
- Boost.Rational;
- Boost.Numeric.Interval;
- platform floating-point rounding control (`fenv`) selected by Boost;
- platform libm behavior used through Boost `rounded_transc_std`;
- explicit pi bracketing constants;
- certificate parser/serializer once external serialization is enabled;
- core finite-arc relation implementation used by the verifier.

Differential tests validate logic above this shared numerical substrate; they do not independently prove the TCB itself. Analytic known-value tests are therefore required for the TCB.
