# Numerical Contract

## Candidate vs certification
Ordinary floating point may propose values but is not proof.

## Certified backend, Phase 1–3
The initial certification substrate is Boost.Numeric.Interval with hardware-directed arithmetic and `rounded_transc_std` for `sin`, `cos`, and `atan`, plus explicit bracketing of rational-to-binary conversion and pi.

This backend and the platform C/libm rounding behavior are part of the Trusted Computing Base (TCB). Certification claims are conditional on those documented TCB semantics.

## Turn representation
A theorem-facing `Turn` is an exact rational multiple of pi. This makes the phase congruences exact and avoids pretending decimal radian approximations satisfy exact modulo-2pi equalities.

## No magic epsilon
Theorem decisions use exact rational comparisons or interval separation. An interval that contains a decision boundary produces `INDETERMINATE`, not a guessed Boolean.

## Certified transcendental policy
All theorem-facing sin/cos/atan evaluations pass through interval operations. `sinc` is represented by a dedicated interval routine; near zero, the routine uses the removable value 1 when the argument is exactly zero and otherwise interval division. The ABCABC theorem domain excludes zero turn itself.
