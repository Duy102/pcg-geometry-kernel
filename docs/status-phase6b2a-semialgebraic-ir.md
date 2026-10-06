# Phase 6B2A — Solver-Neutral Semialgebraic IR

## Scope

Phase 6B2A does **not** claim that a CAD / real-quantifier-elimination
backend has been implemented.

It freezes the theorem-to-backend boundary so that an exact backend does not
have to reinterpret the PCG geometry informally.

## Theorem source

`PCG_General_Fixed_Turn_Intersection_Feasibility_v0.1.tex`

SHA-256:

`625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208`

The theorem gives the complete decision sentence

`exists (p,c) [Network && Positive && Normalize && DistinctVertices && AND_{e<f} not exists q I_ef]`.

For rational multiples of pi, all fixed trigonometric coefficients are
algebraic. Real-closed-field quantifier elimination therefore gives an exact
decision procedure, with CAD as one possible implementation route.

## IR contract

Schema: `pcg-fixed-turn-semialgebraic-ir`

Version: `1.0`

Compiler contract: `PCG-FT-SA-IR-001`

The IR records:

- the exact rational-pi turn data and canonical input SHA-256;
- one edge record for every source edge;
- exact tangent and chord phases;
- turn sign and the theorem's minor / semicircle / major membership branch;
- all quotient-vertex distinctness pairs;
- all unordered finite-arc pairs;
- the allowed common quotient vertices for each pair;
- all fixed prescribed-transversality checks induced by repeated quotient
  symbols;
- outer and pair-local real-variable counts.

The fixed compiler contract defines how each edge record lowers to the theorem
objects `t_e, u_e, lambda_e, o_e, F_e, A_e, B_e` and how every pair record
lowers to `not exists q I_ef`.

## Trust boundary

The kernel can deterministically recompile the IR from its bound input and
reject any mutation of:

- theorem/source identity;
- input binding;
- edge phases or arc class;
- distinctness constraints;
- forbidden pair enumeration;
- allowed common vertices;
- transversality metadata;
- quantifier counts.

The IR has its own SHA-256 digest.

A future Phase 6B2B backend must consume this exact contract. Backend output
must not become a production `REALIZABLE` or `NOT_REALIZABLE` decision
until its evidence has an independently checkable path.

## Why this phase exists

Without a frozen IR, integrating CAD directly into the kernel would mix three
different correctness questions:

1. did PCG compile the theorem correctly?
2. did the backend solve the compiled formula correctly?
3. did the kernel bind the answer to the exact input/theorem version?

Phase 6B2A isolates question 1 and the binding part of question 3 before any
solver is trusted.
