#include <pcg/seifert_flow.hpp>

#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int failures=0;

void check(bool cond, const std::string& name) {
    if (!cond) {
        std::cerr << "FAIL: " << name << "\n";
        ++failures;
    }
}

bool direct_cut_feasible(const pcg::SeifertTurnFlowInput& in) {
    const std::size_t n=in.vertices.size();
    const std::uint64_t total=std::uint64_t{1} << n;
    for (std::uint64_t mask=1;mask<total;++mask) {
        pcg::BigInt P=0,N=0,eps=0,d=0;
        for (std::size_t v=0;v<n;++v) {
            if ((mask>>v)&1U) {
                P+=in.vertices[v].positive_pieces;
                N+=in.vertices[v].negative_pieces;
                eps+=in.vertices[v].orientation;
            }
        }
        for (const auto& e : in.edges) {
            if (e.u==e.v) continue;
            const bool a=((mask>>e.u)&1U)!=0;
            const bool b=((mask>>e.v)&1U)!=0;
            if (a!=b) d+=1;
        }
        if (!(-pcg::BigInt{2}*N-d < pcg::BigInt{2}*eps &&
              pcg::BigInt{2}*eps < pcg::BigInt{2}*P+d))
            return false;
    }
    return true;
}

std::uint64_t rng_state=0x9e3779b97f4a7c15ULL;
std::uint64_t next_u64() {
    rng_state ^= rng_state << 7;
    rng_state ^= rng_state >> 9;
    rng_state ^= rng_state << 8;
    return rng_state;
}

pcg::SeifertTurnFlowInput random_valid_model() {
    pcg::SeifertTurnFlowInput in;
    const std::size_t n=2+(next_u64()%5);
    in.vertices.resize(n);

    for (auto& v : in.vertices) {
        std::uint64_t P=next_u64()%3;
        std::uint64_t N=next_u64()%3;
        if (P+N==0) P=1;
        v.positive_pieces=pcg::BigInt{P};
        v.negative_pieces=pcg::BigInt{N};
        v.orientation=(next_u64()&1U) ? +1 : -1;
    }

    for (std::size_t v=1;v<n;++v)
        in.edges.push_back({v-1,v});

    const std::size_t extras=next_u64()%(n+4);
    for (std::size_t k=0;k<extras;++k) {
        const std::size_t u=next_u64()%n;
        const std::size_t v=(next_u64()%10==0) ? u : next_u64()%n;
        in.edges.push_back({u,v});
    }
    return in;
}

} // namespace

int main() {
    using namespace pcg;

    SeifertTurnFlowInput single_bad{{SeifertVertexCapacity{1,0,+1}},{}};
    auto sb=solve_seifert_turn_flow(single_bad);
    check(sb.status==SeifertFlowStatus::Infeasible,
          "single vertex detects strict positive cut obstruction");
    check(sb.proof==SeifertFlowProof::PositiveCutObstruction,
          "single vertex returns positive cut certificate");
    check(verify_seifert_turn_flow_certificate(sb.certificate).status==SeifertFlowStatus::Infeasible,
          "single vertex obstruction certificate verifies");

    SeifertTurnFlowInput single_good{{SeifertVertexCapacity{2,0,+1}},{}};
    auto sg=solve_seifert_turn_flow(single_good);
    check(sg.status==SeifertFlowStatus::Feasible,
          "single vertex interior turn is feasible");
    check(sg.assurance==ArithmeticAssurance::Exact,
          "Seifert flow witness is exact");
    check(verify_seifert_turn_flow_certificate(sg.certificate).status==SeifertFlowStatus::Feasible,
          "single vertex constructive witness verifies");

    SeifertTurnFlowInput two{{
        SeifertVertexCapacity{1,1,+1},
        SeifertVertexCapacity{1,1,-1}
    },{{0,1}}};
    auto tw=solve_seifert_turn_flow(two);
    check(tw.status==SeifertFlowStatus::Feasible,
          "two-vertex bounded-divergence instance is feasible");
    check(verify_seifert_turn_flow_certificate(tw.certificate).status==SeifertFlowStatus::Feasible,
          "two-vertex exact witness verifies");

    auto two_reversed=two;
    two_reversed.edges[0]={1,0};
    check(canonicalize_seifert_turn_flow_input(two)==
          canonicalize_seifert_turn_flow_input(two_reversed),
          "canonical input is independent of arbitrary edge orientation");
    check(solve_seifert_turn_flow(two_reversed).status==SeifertFlowStatus::Feasible,
          "edge orientation reversal preserves feasibility");

    auto looped=single_good;
    looped.edges={{0,0},{0,0}};
    check(solve_seifert_turn_flow(looped).status==SeifertFlowStatus::Feasible,
          "self-loops do not affect relaxed turn-flow feasibility");

    SeifertTurnFlowInput leaf{{
        SeifertVertexCapacity{0,1,+1},
        SeifertVertexCapacity{3,1,+1},
        SeifertVertexCapacity{2,1,+1}
    },{{0,1},{0,1},{1,2},{1,2}}};
    auto lf=solve_seifert_turn_flow(leaf);
    check(lf.status==SeifertFlowStatus::Infeasible,
          "degree-two all-negative leaf reproduces v0.2 cut obstruction");
    check(lf.proof==SeifertFlowProof::PositiveCutObstruction,
          "leaf obstruction is positive cut budget");
    check(verify_seifert_turn_flow_certificate(lf.certificate).status==SeifertFlowStatus::Infeasible,
          "leaf cut certificate independently verifies");

    auto tampered=tw.certificate;
    tampered.scaled_source_turn[0]+=1;
    auto bad_verify=verify_seifert_turn_flow_certificate(tampered);
    check(bad_verify.status==SeifertFlowStatus::Indeterminate &&
          bad_verify.termination==TerminationReason::BackendFailure,
          "verifier rejects tampered exact flow witness");

    auto wrong_digest=tw.certificate;
    wrong_digest.canonical_input_digest=std::string(64,'0');
    check(verify_seifert_turn_flow_certificate(wrong_digest).status==SeifertFlowStatus::Indeterminate,
          "verifier rejects tampered input digest");

    auto bad_cut=lf.certificate;
    bad_cut.cut_vertices.assign(leaf.vertices.size(),false);
    check(verify_seifert_turn_flow_certificate(bad_cut).status==SeifertFlowStatus::Indeterminate,
          "verifier rejects empty obstruction cut");

    int disagreements=0;
    for (int i=0;i<6000;++i) {
        auto in=random_valid_model();
        const bool expected=direct_cut_feasible(in);
        const auto got=solve_seifert_turn_flow(in);
        const bool actual=got.status==SeifertFlowStatus::Feasible;
        if (got.status==SeifertFlowStatus::Indeterminate ||
            actual!=expected ||
            verify_seifert_turn_flow_certificate(got.certificate).status!=got.status)
            ++disagreements;
    }
    check(disagreements==0,
          "6000/6000 deterministic valid-model cases agree with exact cut oracle");

    if (failures) {
        std::cerr << failures << " Seifert-flow test(s) failed\n";
        return 1;
    }

    std::cout << "Seifert turn-flow tests passed, including 6000-case exact cut cross-check\n";
    return 0;
}
