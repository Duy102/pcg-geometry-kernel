#pragma once

#include <pcg/core.hpp>
#include <cstddef>
#include <string>
#include <vector>

namespace pcg {

struct SeifertVertexCapacity {
    BigInt positive_pieces{0};
    BigInt negative_pieces{0};
    int orientation{+1};
};

struct SeifertEdge {
    std::size_t u{};
    std::size_t v{};
};

struct SeifertTurnFlowInput {
    std::vector<SeifertVertexCapacity> vertices;
    std::vector<SeifertEdge> edges;
};

enum class SeifertFlowStatus {
    Feasible,
    Infeasible,
    Indeterminate
};

enum class SeifertFlowProof {
    InteriorFlowWitness,
    PositiveCutObstruction,
    NegativeCutObstruction,
    None
};

struct SeifertFlowCertificate {
    explicit SeifertFlowCertificate(const SeifertTurnFlowInput& in) : input(in) {}

    std::string schema{"pcg-seifert-turn-flow-certificate"};
    std::string schema_version{"1.0"};
    std::string theorem_id{"PCG-SEIFERT-TURN-FLOW-THM-001"};
    SourceDigest source_digest{
        "SHA-256",
        "1a60d9ea904775dfd59c345915e1cb00fb67f13dbc46a029a9a7c8846428daa5"
    };
    std::string canonical_input_digest;

    SeifertTurnFlowInput input;
    SeifertFlowStatus status{SeifertFlowStatus::Indeterminate};
    SeifertFlowProof proof{SeifertFlowProof::None};

    BigInt denominator{0};
    std::vector<BigInt> scaled_edge_flow;
    std::vector<BigInt> scaled_source_turn;
    std::vector<bool> cut_vertices;
};

struct SeifertFlowResult {
    explicit SeifertFlowResult(const SeifertFlowCertificate& c) : certificate(c) {}

    SeifertFlowStatus status{SeifertFlowStatus::Indeterminate};
    SeifertFlowProof proof{SeifertFlowProof::None};
    ArithmeticAssurance assurance{ArithmeticAssurance::None};
    TerminationReason termination{TerminationReason::BackendFailure};
    SeifertFlowCertificate certificate;
};

std::string canonicalize_seifert_turn_flow_input(const SeifertTurnFlowInput& input);
std::string serialize_seifert_turn_flow_certificate(const SeifertFlowCertificate& certificate);

SeifertFlowResult solve_seifert_turn_flow(const SeifertTurnFlowInput& input);
SeifertFlowResult verify_seifert_turn_flow_certificate(const SeifertFlowCertificate& certificate);

} // namespace pcg
