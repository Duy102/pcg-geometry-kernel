# PCG Geometry Kernel

**Status: Phase 0 audit checkpoint. No production solver is implemented yet.**

This public repository holds the PCG audit package approved for publication.
The immediate build target is a certified ABCABC vertical slice (Phases 0–3),
with an independently checked certificate verifier.

## Audit package

[Download PCG_Geometry_Kernel_Phase0_Audit.zip](PCG_Geometry_Kernel_Phase0_Audit.zip)

Archive SHA-256:

```text
2633531cc90a083073fcd4c21693d535086d275d27266e6ac9d16f53296a6880
```

The archive contains 43 files under `pcg-geometry-kernel/`:

- The supplied production charter and original source snapshots.
- Draft product, numerical, degeneracy, canonicalization, provenance and TCB contracts.
- A source inventory and 12 stable statement records with SHA-256 digests.
- A blocker register, an audit report and the Phase 0–3 implementation plan.
- Five source-witness candidates and a C++ backend readiness probe.

Extract the archive, then start with `docs/audit-report.md` and
`docs/blockers.md`. Source and statement hashes identify the exact source
versions; they do not establish theorem truth.

## Verified at this checkpoint

- 23 source snapshots: byte digests and sizes checked.
- 12 statement records: unique identifiers and statement digests checked.
- ABCABC pair partition: exactly 15 = 3 + 6 + 6.
- Five source-witness candidates: exact turn-range and phase-congruence checks only.

No complete geometric witness, numerical certification, golden corpus, certificate
verifier, production solver or GitHub Actions result is claimed.

## Open blockers

1. The original execution environment lacks the MPFR/GMP development files;
   the readiness probe failed at missing `mpfr.h`.
2. The precise executable interpretation of ABCABC genericity must be frozen.
3. Global sharp-refinement source status differs across supplied manuscripts;
   those claims are quarantined pending model/source reconciliation. This is
   not a demonstrated counterexample to the separate fixed-turn ABCABC theorem.

The implementation language selected is C++20; the intended numerical stack is
GMP/MPFR with rigorously validated enclosures. The required backend is designed,
not yet implemented and tested.

## Publication record

The public repository and this archive were published after the audit was
packaged on 2026-10-06. The archive's statements that no repository had been
created describe the earlier audit checkpoint; they are retained as historical
records. Its mathematical and build blockers remain open.

Public visibility does not grant a software or manuscript license.
No license has been added.
