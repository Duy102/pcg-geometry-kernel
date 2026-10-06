#pragma once

#include <pcg/core.hpp>
#include <cstddef>
#include <string>
#include <vector>

namespace pcg {

struct NetworkClosureInput {
    std::size_t vertex_count{};
    std::vector<std::size_t> trace_vertices;
    std::vector<Turn> turns;
};

struct NetworkCycleBasis {
    std::vector<std::vector<int>> rows;
    std::vector<PiRational> chord_phase_pi;
};

enum class NetworkClosureStatus {
    Closed,
    NotClosed,
    Indeterminate
};

enum class NetworkClosureProof {
    ExactPositiveRigidKernel,
    ExactFullRankObstruction,
    ExactRigidSignObstruction,
    HigherDimensionalKernel,
    CyclotomicOrderLimit,
    SignCertificationLimit,
    BackendFailure,
    None
};

struct AlgebraicChordValue {
    // Coefficients in the power basis 1,zeta,...,zeta^(phi(n)-1)
    // of Q(zeta_n), where n = cyclotomic_order.
    std::vector<Rational> coefficients;
};

struct NetworkClosureCertificate {
    explicit NetworkClosureCertificate(const NetworkClosureInput& in) : input(in) {}

    std::string schema{"pcg-network-closure-certificate"};
    std::string schema_version{"1.0"};
    std::string theorem_id{"PCG-NETWORK-CLOSURE-THM-001"};
    SourceDigest source_digest{
        "SHA-256",
        "6d1f8fc39b3bb34385956c4d2c0e77e222260f60d2f59fde860d6c1964aa8546"
    };
    std::string canonical_input_digest;

    NetworkClosureInput input;
    NetworkClosureStatus status{NetworkClosureStatus::Indeterminate};
    NetworkClosureProof proof{NetworkClosureProof::None};

    unsigned cyclotomic_order{};
    std::size_t exact_rank{};
    std::vector<AlgebraicChordValue> rigid_kernel;
};

struct NetworkClosureResult {
    explicit NetworkClosureResult(const NetworkClosureCertificate& c) : certificate(c) {}

    NetworkClosureStatus status{NetworkClosureStatus::Indeterminate};
    NetworkClosureProof proof{NetworkClosureProof::None};
    ArithmeticAssurance assurance{ArithmeticAssurance::None};
    TerminationReason termination{TerminationReason::BackendFailure};
    NetworkClosureCertificate certificate;
};

NetworkCycleBasis build_network_cycle_basis(const NetworkClosureInput& input);
std::string canonicalize_network_closure_input(const NetworkClosureInput& input);
std::string serialize_network_closure_certificate(const NetworkClosureCertificate& certificate);

NetworkClosureResult solve_network_closure(const NetworkClosureInput& input);
NetworkClosureResult verify_network_closure_certificate(const NetworkClosureCertificate& certificate);

} // namespace pcg
