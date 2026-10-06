# Phase 4 Status — Seifert Turn Flow

Status: production-candidate module on branch phase4-seifert-turn-flow.

Implemented:

- exact arbitrary-precision input capacities;
- fixed Seifert orientation vector;
- multigraph support with loops ignored in net flow as required by the source;
- exact common interior shrink;
- lower-bound circulation reduction;
- Dinic max-flow over arbitrary-precision integer capacities;
- constructive exact FEASIBLE witness;
- exact cut INFEASIBLE witness;
- certificate binding to canonical input and theorem source SHA-256;
- independent algebraic certificate verification without rerunning max-flow;
- 6000-case deterministic exact cut cross-check;
- GCC/Clang Debug and Release plus ASan/UBSan CI.

Scientific boundary:

Phase 4 decides only the relaxed Seifert turn-flow layer. Network Closure and trace-faithful Euclidean realization are not part of this module.

Source-text note:

PCG_Seifert_Turn_Flow_v0.2.tex states the strict cut family for all subsets U, but the empty subset would make the displayed strict inequality impossible. The implementation uses nonempty U, consistent with the circulation proof, and records this discrepancy in docs/seifert-turn-flow-contract.md.