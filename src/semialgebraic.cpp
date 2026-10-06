#include <pcg/semialgebraic.hpp>

#include <algorithm>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pcg {
namespace {

constexpr const char* kSchema="pcg-fixed-turn-semialgebraic-ir";
constexpr const char* kSchemaVersion="1.0";
constexpr const char* kCompilerContract="PCG-FT-SA-IR-001";
constexpr const char* kTheoremId="PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001";
constexpr const char* kSourceDigest=
    "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208";

Rational abs_rational(Rational x) {
    return x < Rational(BigInt{0}) ? -x : x;
}

SemialgebraicArcClass classify_arc(const Turn& t) {
    const Rational a=abs_rational(t.pi.value);
    const Rational one(BigInt{1});
    if (a < one) return SemialgebraicArcClass::Minor;
    if (a == one) return SemialgebraicArcClass::Semicircle;
    return SemialgebraicArcClass::Major;
}

int turn_sign(const Turn& t) {
    return t.pi.value < Rational(BigInt{0}) ? -1 : 1;
}

std::vector<std::size_t> shared_endpoint_vertices(
    const NetworkClosureInput& in,std::size_t e,std::size_t f) {
    const std::size_t m=in.trace_vertices.size();
    const std::size_t et=in.trace_vertices[e];
    const std::size_t eh=in.trace_vertices[(e+1)%m];
    const std::size_t ft=in.trace_vertices[f];
    const std::size_t fh=in.trace_vertices[(f+1)%m];

    std::set<std::size_t> a{et,eh};
    std::set<std::size_t> b{ft,fh};
    std::vector<std::size_t> out;
    std::set_intersection(a.begin(),a.end(),b.begin(),b.end(),
                          std::back_inserter(out));
    return out;
}

void append_size_vector(std::ostringstream& os,
                        const std::vector<std::size_t>& xs) {
    os << "[";
    for (std::size_t i=0;i<xs.size();++i) {
        if (i) os << ",";
        os << xs[i];
    }
    os << "]";
}

void append_input(std::ostringstream& os,const NetworkClosureInput& in) {
    os << "{\"vertex_count\":" << in.vertex_count
       << ",\"trace_vertices\":";
    append_size_vector(os,in.trace_vertices);
    os << ",\"turns_pi\":[";
    for (std::size_t i=0;i<in.turns.size();++i) {
        if (i) os << ",";
        os << "\"" << rational_string(in.turns[i].pi) << "\"";
    }
    os << "]}";
}

bool metadata_ok(const FixedTurnSemialgebraicProblem& p) {
    return p.schema==kSchema &&
           p.schema_version==kSchemaVersion &&
           p.compiler_contract==kCompilerContract &&
           p.theorem_id==kTheoremId &&
           p.source_digest.algorithm=="SHA-256" &&
           p.source_digest.value==kSourceDigest;
}

} // namespace

FixedTurnSemialgebraicProblem compile_fixed_turn_semialgebraic_problem(
    const NetworkClosureInput& input) {
    // Reuse the exact Network Closure front-end for graph/input validation and
    // exact propagated chord phases.
    const auto basis=build_network_cycle_basis(input);

    FixedTurnSemialgebraicProblem out;
    out.input=input;
    out.canonical_input_digest=sha256_hex(
        canonicalize_network_closure_input(input));
    out.root_vertex=0;

    const std::size_t m=input.turns.size();
    const std::size_t n=input.vertex_count;

    out.edges.reserve(m);
    PiRational tangent{0,1};
    for (std::size_t e=0;e<m;++e) {
        const auto expected_chord=tangent+input.turns[e].pi/2;
        if (!(expected_chord==basis.chord_phase_pi[e]))
            throw std::runtime_error("semialgebraic compiler phase mismatch");

        out.edges.push_back(SemialgebraicEdgeSpec{
            e,
            input.trace_vertices[e],
            input.trace_vertices[(e+1)%m],
            tangent,
            expected_chord,
            input.turns[e].pi,
            turn_sign(input.turns[e]),
            classify_arc(input.turns[e])
        });
        tangent=tangent+input.turns[e].pi;
    }

    for (std::size_t v=0;v<n;++v)
        for (std::size_t w=v+1;w<n;++w)
            out.distinct_vertex_pairs.push_back({v,w});

    for (std::size_t e=0;e<m;++e) {
        for (std::size_t f=e+1;f<m;++f) {
            out.forbidden_pairs.push_back(
                SemialgebraicPairSpec{
                    e,f,shared_endpoint_vertices(input,e,f)});
        }
    }

    std::vector<std::vector<std::size_t>> occurrences(n);
    for (std::size_t e=0;e<m;++e)
        occurrences[input.trace_vertices[e]].push_back(e);

    for (std::size_t v=0;v<n;++v) {
        const auto& occ=occurrences[v];
        for (std::size_t i=0;i<occ.size();++i) {
            for (std::size_t j=i+1;j<occ.size();++j) {
                const auto delta=
                    out.edges[occ[j]].tangent_phase_pi-
                    out.edges[occ[i]].tangent_phase_pi;
                out.prescribed_transversality.push_back(
                    SemialgebraicTransversalitySpec{
                        v,occ[i],occ[j],delta,sin_pi_sign(delta)!=0});
            }
        }
    }

    out.outer_real_variable_count=2*n+m;
    out.pair_local_real_variable_count=2*out.forbidden_pairs.size();
    return out;
}

std::string serialize_fixed_turn_semialgebraic_problem(
    const FixedTurnSemialgebraicProblem& p) {
    std::ostringstream os;
    os << "{\"schema\":\"" << p.schema
       << "\",\"schema_version\":\"" << p.schema_version
       << "\",\"compiler_contract\":\"" << p.compiler_contract
       << "\",\"theorem_id\":\"" << p.theorem_id
       << "\",\"source_digest\":{\"algorithm\":\"" << p.source_digest.algorithm
       << "\",\"value\":\"" << p.source_digest.value
       << "\"},\"canonical_input_digest\":\"" << p.canonical_input_digest
       << "\",\"input\":";
    append_input(os,p.input);
    os << ",\"root_vertex\":" << p.root_vertex;

    os << ",\"edges\":[";
    for (std::size_t i=0;i<p.edges.size();++i) {
        if (i) os << ",";
        const auto& e=p.edges[i];
        os << "{\"edge\":" << e.edge
           << ",\"tail\":" << e.tail
           << ",\"head\":" << e.head
           << ",\"tangent_phase_pi\":\"" << rational_string(e.tangent_phase_pi)
           << "\",\"chord_phase_pi\":\"" << rational_string(e.chord_phase_pi)
           << "\",\"turn_pi\":\"" << rational_string(e.turn_pi)
           << "\",\"turn_sign\":" << e.turn_sign
           << ",\"arc_class\":" << static_cast<int>(e.arc_class)
           << "}";
    }
    os << "]";

    os << ",\"distinct_vertex_pairs\":[";
    for (std::size_t i=0;i<p.distinct_vertex_pairs.size();++i) {
        if (i) os << ",";
        os << "[" << p.distinct_vertex_pairs[i].first
           << "," << p.distinct_vertex_pairs[i].second << "]";
    }
    os << "]";

    os << ",\"forbidden_pairs\":[";
    for (std::size_t i=0;i<p.forbidden_pairs.size();++i) {
        if (i) os << ",";
        const auto& x=p.forbidden_pairs[i];
        os << "{\"e\":" << x.e << ",\"f\":" << x.f
           << ",\"allowed_common_vertices\":";
        append_size_vector(os,x.allowed_common_vertices);
        os << "}";
    }
    os << "]";

    os << ",\"prescribed_transversality\":[";
    for (std::size_t i=0;i<p.prescribed_transversality.size();++i) {
        if (i) os << ",";
        const auto& x=p.prescribed_transversality[i];
        os << "{\"vertex\":" << x.vertex
           << ",\"occurrence_a\":" << x.occurrence_a
           << ",\"occurrence_b\":" << x.occurrence_b
           << ",\"tangent_delta_pi\":\"" << rational_string(x.tangent_delta_pi)
           << "\",\"transverse\":" << (x.transverse?"true":"false")
           << "}";
    }
    os << "]";

    os << ",\"outer_real_variable_count\":" << p.outer_real_variable_count
       << ",\"pair_local_real_variable_count\":" << p.pair_local_real_variable_count
       << ",\"quantified_formula\":\"exists(p,c)[Network&&Positive&&Normalize&&DistinctVertices&&AND_pair(not exists(q) I_ef)]\""
       << "}";
    return os.str();
}

std::string fixed_turn_semialgebraic_problem_digest(
    const FixedTurnSemialgebraicProblem& problem) {
    return sha256_hex(serialize_fixed_turn_semialgebraic_problem(problem));
}

SemialgebraicProblemValidation verify_fixed_turn_semialgebraic_problem(
    const FixedTurnSemialgebraicProblem& problem) {
    try {
        if (!metadata_ok(problem))
            return {SemialgebraicProblemValidationStatus::InvalidMetadata,
                    "schema, compiler contract, theorem id, or source digest mismatch"};

        const auto expected_input_digest=sha256_hex(
            canonicalize_network_closure_input(problem.input));
        if (problem.canonical_input_digest!=expected_input_digest)
            return {SemialgebraicProblemValidationStatus::InvalidInputBinding,
                    "canonical input digest mismatch"};

        const auto expected=compile_fixed_turn_semialgebraic_problem(problem.input);
        if (serialize_fixed_turn_semialgebraic_problem(problem)!=
            serialize_fixed_turn_semialgebraic_problem(expected))
            return {SemialgebraicProblemValidationStatus::InvalidStructure,
                    "IR differs from deterministic theorem compiler output"};

        return {SemialgebraicProblemValidationStatus::Valid,""};
    } catch (const std::exception& e) {
        return {SemialgebraicProblemValidationStatus::InvalidStructure,e.what()};
    } catch (...) {
        return {SemialgebraicProblemValidationStatus::InvalidStructure,
                "unknown semialgebraic IR verification failure"};
    }
}

} // namespace pcg
