#pragma once

#include <pcg/network_closure.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace pcg {

enum class SemialgebraicArcClass {
    Minor,
    Semicircle,
    Major
};

struct SemialgebraicEdgeSpec {
    std::size_t edge{};
    std::size_t tail{};
    std::size_t head{};
    PiRational tangent_phase_pi;
    PiRational chord_phase_pi;
    PiRational turn_pi;
    int turn_sign{};
    SemialgebraicArcClass arc_class{SemialgebraicArcClass::Minor};
};

struct SemialgebraicPairSpec {
    std::size_t e{};
    std::size_t f{};
    std::vector<std::size_t> allowed_common_vertices;
};

struct SemialgebraicTransversalitySpec {
    std::size_t vertex{};
    std::size_t occurrence_a{};
    std::size_t occurrence_b{};
    PiRational tangent_delta_pi;
    bool transverse{};
};

struct FixedTurnSemialgebraicProblem {
    std::string schema{"pcg-fixed-turn-semialgebraic-ir"};
    std::string schema_version{"1.0"};
    std::string compiler_contract{"PCG-FT-SA-IR-001"};
    std::string theorem_id{"PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001"};
    SourceDigest source_digest{
        "SHA-256",
        "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208"
    };
    std::string canonical_input_digest;

    NetworkClosureInput input;
    std::size_t root_vertex{};

    // The contract semantics are:
    //   exists p,c [
    //     Network(p,c) && Positive(c) && Normalize(c)
    //     && DistinctVertices(p)
    //     && AND_{e<f} ! exists q I_ef(p,c,q)
    //   ].
    //
    // Edge specs carry the exact rational-pi data needed to derive
    // t_e, u_e, lambda_e, o_e, F_e, A_e and B_e.
    std::vector<SemialgebraicEdgeSpec> edges;
    std::vector<std::pair<std::size_t,std::size_t>> distinct_vertex_pairs;
    std::vector<SemialgebraicPairSpec> forbidden_pairs;
    std::vector<SemialgebraicTransversalitySpec> prescribed_transversality;

    std::size_t outer_real_variable_count{};
    std::size_t pair_local_real_variable_count{};
};

enum class SemialgebraicProblemValidationStatus {
    Valid,
    InvalidMetadata,
    InvalidInputBinding,
    InvalidStructure
};

struct SemialgebraicProblemValidation {
    SemialgebraicProblemValidationStatus status{
        SemialgebraicProblemValidationStatus::InvalidStructure};
    std::string detail;
};

FixedTurnSemialgebraicProblem compile_fixed_turn_semialgebraic_problem(
    const NetworkClosureInput& input);

std::string serialize_fixed_turn_semialgebraic_problem(
    const FixedTurnSemialgebraicProblem& problem);

std::string fixed_turn_semialgebraic_problem_digest(
    const FixedTurnSemialgebraicProblem& problem);

SemialgebraicProblemValidation verify_fixed_turn_semialgebraic_problem(
    const FixedTurnSemialgebraicProblem& problem);

} // namespace pcg
