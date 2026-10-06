# Phase 4 — General Network Positive-Kernel Closure

This phase extends the existing Network Closure implementation beyond the projectively rigid (one-dimensional kernel) case.

## Theorem basis

Source: `PCG_Network_Closure_Theorem_v0.1.tex`

SHA-256:

`6d1f8fc39b3bb34385956c4d2c0e77e222260f60d2f59fde860d6c1964aa8546`

The theorem requires deciding whether the exact direction-network matrix `N` has a strictly positive kernel vector:

```text
exists c > 0 such that N c = 0
```

## Production extension

For exact cyclotomic systems within the existing order cap:

- rank = number of columns: exact full-rank obstruction;
- kernel dimension = 1: retain the existing rigid-kernel path;
- kernel dimension > 1: construct an exact kernel basis and solve the homogeneous positivity problem by Fourier–Motzkin elimination over the real cyclotomic field.

Strict positivity is converted using homogeneity:

```text
K z > 0  iff  there exists a scaling with K z >= 1.
```

Elimination arithmetic is exact in the cyclotomic field. Ordering decisions use the existing certified real-sign oracle. Therefore uncertain sign separation remains `INDETERMINATE`; it is never guessed.

To prevent uncontrolled exponential elimination growth, the implementation has a finite inequality-count cap. Hitting that cap returns `INDETERMINATE` with `EliminationComplexityLimit`.

## Certificates

A feasible higher-dimensional instance stores an explicit algebraic positive kernel witness and the verifier checks:

- theorem/source/input binding;
- exact matrix-kernel equality;
- strict positivity of every chord magnitude.

An infeasible higher-dimensional instance uses `CertifiedGeneralConeObstruction`; the verifier independently reruns the deterministic exact/certified elimination path.

## Regression cases

- 4-cycle with four turns `pi/2`: rank 2, kernel dimension 2, positive closure exists.
- 4-cycle with four turns `pi/4`: rank 2, kernel dimension 2, all directions lie in one open half-plane, so positive closure is impossible.
