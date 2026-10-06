#pragma once

#include <pcg/geometry.hpp>
#include <pcg/network_closure.hpp>
#include <cstddef>
#include <limits>
#include <string>

namespace pcg {

enum class RigidPostClosureProof {
    AllPairsClear,
    NetworkClosureObstruction,
    ExactVertexCollision,
    CertifiedUnintendedIntersection,
    NonTransversePrescribedPassage,
    UnsupportedTraceMultiplicity,
    HigherDimensionalClosure,
    PairCertificationLimit,
    BackendFailure,
    None
};

struct RigidPostClosureCertificate {
    explicit RigidPostClosureCertificate(const NetworkClosureCertificate& c)
        : closure_certificate(c) {}

    std::string schema{"pcg-rigid-post-closure-certificate"};
    std::string schema_version{"1.0"};
    std::string theorem_id{"PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001"};
    SourceDigest source_digest{
        "SHA-256",
        "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208"
    };
    std::string closure_binding_digest;

    NetworkClosureCertificate closure_certificate;
    Decision decision{Decision::Indeterminate};
    RigidPostClosureProof proof{RigidPostClosureProof::None};
    ArithmeticAssurance assurance{ArithmeticAssurance::None};
    TerminationReason termination{TerminationReason::BackendFailure};

    std::size_t witness_a{std::numeric_limits<std::size_t>::max()};
    std::size_t witness_b{std::numeric_limits<std::size_t>::max()};
    std::size_t checked_pairs{};
};

struct RigidPostClosureResult {
    explicit RigidPostClosureResult(const RigidPostClosureCertificate& c)
        : certificate(c) {}

    Decision decision{Decision::Indeterminate};
    RigidPostClosureProof proof{RigidPostClosureProof::None};
    ArithmeticAssurance assurance{ArithmeticAssurance::None};
    TerminationReason termination{TerminationReason::BackendFailure};
    RigidPostClosureCertificate certificate;
};

std::string serialize_rigid_post_closure_certificate(
    const RigidPostClosureCertificate& certificate);

RigidPostClosureResult solve_projectively_rigid_post_closure(
    const NetworkClosureCertificate& closure_certificate);

RigidPostClosureResult verify_projectively_rigid_post_closure_certificate(
    const RigidPostClosureCertificate& certificate);

} // namespace pcg
