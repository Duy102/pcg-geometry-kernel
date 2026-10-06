# Phase 6B2C — Exact QEPCAD / RCF Backend Adapter

## Scope

Phase 6B2C converts the verified Phase 6B2B exact-algebraic polynomial AST into
a rational-coefficient first-order formula suitable for an exact
real-closed-field (RCF) quantifier-elimination backend.

The first concrete interchange target is **QEPCAD B**.

This phase is an adapter and trust-boundary phase. It does **not** execute
QEPCAD and it does **not** promote an external SAT / UNSAT answer into a PCG
production decision.

## Why another lowering is needed

Phase 6B2B intentionally keeps theorem coefficients in exact symbolic form:

- rational numbers;
- `sin(pi*r)`;
- `cos(pi*r)`;
- addition;
- multiplication;
- negation;
- inversion.

A CAD backend expects polynomial formulas over real variables with rational
coefficients. Phase 6B2C therefore replaces every algebraic coefficient by
real variables plus exact polynomial definitions.

No theorem coefficient is converted to ordinary floating point.

## Exact rational-pi trigonometric constants

For a reduced rational phase

`r = n/d`

and

`c = cos(pi*r)`

the adapter uses the exact Chebyshev relation

`T_d(c) - (-1)^n = 0`.

That equation alone can have several distinct real roots, so the adapter also
selects the intended root exactly.

It computes the target root rank from `n mod 2d`, then constructs a rational
isolating interval using exact Sturm root counting. The interval is accepted
only when it contains exactly one distinct root of the Chebyshev relation.

This isolation path uses arbitrary-precision rational arithmetic only; it does
not depend on `libm`, decimal approximations, or a tolerance epsilon.

For

`s = sin(pi*r)`

the adapter adds

`c^2 + s^2 - 1 = 0`

and the exact sign of `s`, obtained from the rational-pi sign logic already
in the PCG core.

Together, the cosine root isolation, unit-circle relation, and exact sine sign
select the intended pair `(c,s)` uniquely.

## Composite algebraic coefficients

Every non-atomic algebraic expression receives a deterministic auxiliary real
variable.

The adapter emits exact defining equations:

- add: `z = a + b + ...`;
- multiply: `z = a b ...`;
- negate: `z = -a`;
- inverse: `z a = 1`.

The inverse equation fails closed if its denominator were zero. For the current
PCG fixed-turn lowering, the generated inverse is based on
`2 sin(tau/2)`, and the theorem input domain `0 < |tau| < 2*pi` excludes
zero.

## Quantifier structure

The algebraic-coefficient variables are existentially quantified outside the
Phase 6B2B theorem sentence.

The Phase 6B2B bad-pair form

`not exists (q_x,q_y) I_ef`

is prenexed for QEPCAD as universal pair-local variables.

Thus a closed PCG decision instance is exported with zero free variables and a
prefix of the form

`(E algebraic constants)(E p,c)(A q_01)(A q_02)...`.

The quantifier-free matrix contains only rational polynomial equalities,
inequalities, Boolean conjunction / disjunction, and negated atomic formulas.

## QEPCAD syntax

The exporter follows QEPCAD B's documented input structure:

1. description in square brackets;
2. ordered variable list;
3. number of free variables;
4. prenex formula terminated by a period;
5. `finish`.

For PCG decision instances the number of free variables is `0`.

QEPCAD's documentation explicitly permits polynomial atomic formulas with
rational coefficients and the Boolean operators `/\\`, `\\/`, and
`~`.

Reference:

`https://www.usna.edu/Users/cs/wcbrown/qepcad/B/user/EnterForm.html`

## Canonical binding

Schema:

`pcg-qepcad-rcf-program`

Version:

`1.0`

Backend contract:

`PCG-FT-QE-QEPCAD-001`

The adapter package binds to:

- the fixed-turn theorem id and source SHA-256;
- the exact Network Closure input;
- the exact SHA-256 of the verified Phase 6B2B program;
- all algebraic variables and operation definitions;
- every Chebyshev polynomial;
- every rational isolating interval;
- the full rational RCF formula.

The verifier recompiles Phase 6B2B and Phase 6B2C deterministically and rejects
metadata, source binding, isolation, degree, atom-count, or structural
mutations.

## Resource guard

The present adapter accepts rational-pi trigonometric phases whose reduced
denominator is at most

`256`.

This is an explicit implementation resource guard, aligned with the existing
exact cyclotomic closure substrate. It is **not** a mathematical limit of the
PCG decidability theorem.

Inputs beyond this adapter limit fail closed instead of silently using
approximate coefficients.

## What Phase 6B2C does not claim

Phase 6B2C does not claim:

- that QEPCAD has been installed in the production CI image;
- that every exported instance is computationally practical for CAD;
- that a QEPCAD textual answer is itself a PCG certificate;
- that the general PCG solver now returns global `REALIZABLE` /
  `NOT_REALIZABLE`.

Those claims require an execution and evidence boundary.

## Next phase

Phase 6B2D should add:

1. a pinned QEPCAD execution environment;
2. deterministic invocation and resource limits;
3. parsing of true / false / timeout / backend failure;
4. independent binding of the solver transcript to the Phase 6B2C digest;
5. a policy that only verified exact backend results can enter the production
   PCG decision path, with all other cases remaining `INDETERMINATE`.
