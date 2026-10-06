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
