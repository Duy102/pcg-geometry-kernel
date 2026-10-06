# Degeneracy Policy

| Case | Phase 0–3 policy |
|---|---|
| `tau = 0` | `DomainError::ZeroTurnOutsideTheoremDomain` |
| `|tau| >= 2*pi` | `DomainError::FullOrOverTurnOutsideTheoremDomain` |
| supporting circles required distinct but certification cannot separate centers | `INDETERMINATE` |
| second circle intersection is tangent/zero within interval uncertainty | `INDETERMINATE` |
| exact positive-closure sine is zero | certified `NOT_REALIZABLE` obstruction |
| finite-arc membership boundary cannot be separated | `INDETERMINATE` |
| unsupported trace family | `UNSUPPORTED` at top-level dispatcher (not used by specialized ABCABC entrypoint) |
