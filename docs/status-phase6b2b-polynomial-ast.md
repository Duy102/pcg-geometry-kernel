# Phase 6B2B — Exact-Algebraic Polynomial AST

## Scope

Phase 6B2B lowers the verified Phase 6B2A solver-neutral IR into an explicit
first-order polynomial formula AST.

It still does **not** embed or trust a CAD / real-quantifier-elimination
backend.

## Theorem fidelity

The lowering implements the exact fixed-turn sentence

`exists (p,c) [Network && Positive && Normalize && DistinctVertices && AND_{e<f} not exists q I_ef]`.

For every edge it constructs:

- exact chord direction `u_e=(cos alpha_e,sin alpha_e)`;
- exact signed-radius coefficient
  `lambda_e = 1/(2 sin(tau_e/2))`;
- supporting-circle center
  `o_e = p_tail + lambda_e c_e J t_e`;
- the quadratic supporting-circle equation `F_e(q)=0`;
- the oriented cross-product polynomials `A_e(q)` and `B_e(q)`;
- the theorem's minor / semicircle / major Boolean membership branch.

For every unordered edge pair it emits:

`not exists (q_x,q_y) I_ef`

and excludes each allowed common quotient vertex with a strict squared-distance
inequality.

## Exact coefficient contract

No theorem coefficient is converted to an ordinary floating-point number.

Coefficients are represented by an exact algebraic expression AST generated
from:

- arbitrary-precision rational numbers;
- `sin(pi*r)` and `cos(pi*r)` for exact rational `r`;
- addition;
- multiplication;
- negation;
- inversion.

Because the input contract has `0 < |tau_e| < 2*pi`, every
`sin(tau_e/2)` denominator used by `lambda_e` is nonzero.

For rational `r`, the sine and cosine atoms denote real algebraic numbers
derived from roots of unity. A later backend adapter may map these atoms into
its native exact real-algebraic-number representation.

## Degree contract

All polynomial atoms are checked during lowering and verification.

The maximum permitted degree is:

`2`.

This covers:

- linear network equations;
- linear positivity / normalization;
- quadratic vertex distinctness;
- quadratic supporting-circle equations;
- quadratic oriented cross-product membership conditions;
- quadratic exclusion of allowed common endpoints.

Any future compiler change that exceeds degree two is rejected by the
Phase 6B2B verifier.

## Trust boundary

Schema:

`pcg-exact-semialgebraic-program`

Version:

`1.0`

Lowering contract:

`PCG-FT-SA-POLY-AST-001`

Each program binds to:

- the theorem/source SHA-256;
- the exact Network Closure input;
- the SHA-256 of its verified Phase 6B2A source IR.

The verifier recompiles Phase 6B2A and Phase 6B2B deterministically and rejects
metadata, binding, formula, variable, atom-count, or degree mutations.

## Remaining backend work

Phase 6B2C should add one backend adapter that maps the exact algebraic
coefficient AST and polynomial formula AST into a real-closed-field decision
engine.

A backend result must not become a production `REALIZABLE` or
`NOT_REALIZABLE` decision until its result/evidence path is independently
checked or otherwise reduced to a trusted exact decision boundary.
