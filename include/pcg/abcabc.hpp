#pragma once
#include <pcg/geometry.hpp>
#include <array>
#include <string>

namespace pcg {

struct ABCABCInput {
    std::array<Turn, 6> turns;
};

enum class ABCABCProofReason {
    AllConditionsSatisfied,
    PhaseCongruenceObstruction,
    PositiveClosureObstruction,
    ExtraHitObstruction,
    NumericalIndeterminacy
};

struct ABCABCCertificate {
    explicit ABCABCCertificate(const ABCABCInput& in) : input(in) {}
    std::string schema{"pcg-abcabc-certificate"};
    std::string schema_version{"1.0"};
    std::string theorem_id{"PCG-ABCABC-THM-001"};
    SourceDigest source_digest{"SHA-256", "2aa832a00031035d964c0d9362e4409861e19b97a9d764f3123a42edef400beb"};
    std::string canonical_input_digest;
    ABCABCInput input;
    ABCABCProofReason reason{ABCABCProofReason::NumericalIndeterminacy};
    ProofKind proof_kind{ProofKind::None};
    ArithmeticAssurance assurance{ArithmeticAssurance::None};
};

struct ABCABCResult {
    explicit ABCABCResult(const ABCABCCertificate& cert) : certificate(cert) {}
    ABCABCResult(Decision d, ProofKind k, ArithmeticAssurance a, TerminationReason t, const ABCABCCertificate& cert)
        : decision(d), proof_kind(k), assurance(a), termination(t), certificate(cert) {}
    Decision decision{Decision::Indeterminate};
    ProofKind proof_kind{ProofKind::None};
    ArithmeticAssurance assurance{ArithmeticAssurance::None};
    TerminationReason termination{TerminationReason::Completed};
    ABCABCCertificate certificate;
};

std::string canonicalize_abcabc_input(const ABCABCInput& input);
std::string serialize_abcabc_certificate(const ABCABCCertificate& certificate);
ABCABCResult solve_abcabc(const ABCABCInput& input);
ABCABCResult verify_abcabc_certificate(const ABCABCCertificate& certificate);

} // namespace pcg
