#include <pcg/seifert_flow.hpp>

#include <algorithm>
#include <climits>
#include <cstdint>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pcg {
namespace {

const char* kSchema = "pcg-seifert-turn-flow-certificate";
const char* kSchemaVersion = "1.0";
const char* kTheoremId = "PCG-SEIFERT-TURN-FLOW-THM-001";
const char* kSourceDigest = "1a60d9ea904775dfd59c345915e1cb00fb67f13dbc46a029a9a7c8846428daa5";

struct DinicEdge {
    int to{};
    std::size_t rev{};
    BigInt cap{0};
};

class Dinic {
public:
    explicit Dinic(std::size_t n) : graph_(n), level_(n), it_(n) {}

    std::size_t add_edge(int from, int to, const BigInt& cap) {
        if (cap < 0) throw std::invalid_argument("negative max-flow capacity");
        const std::size_t fi = graph_[from].size();
        const std::size_t ri = graph_[to].size();
        graph_[from].push_back(DinicEdge{to,ri,cap});
        graph_[to].push_back(DinicEdge{from,fi,BigInt{0}});
        return fi;
    }

    BigInt max_flow(int s, int t, const BigInt& limit) {
        BigInt flow = 0;
        while (flow < limit && bfs(s,t)) {
            std::fill(it_.begin(),it_.end(),0);
            while (flow < limit) {
                const BigInt remaining = limit-flow;
                BigInt pushed=dfs(s,t,remaining);
                if (pushed == 0) break;
                flow += pushed;
            }
        }
        return flow;
    }

    std::vector<bool> reachable(int s) const {
        std::vector<bool> seen(graph_.size(),false);
        std::queue<int> q;
        seen[s]=true;
        q.push(s);
        while(!q.empty()) {
            const int v=q.front();
            q.pop();
            for (const auto& e : graph_[v]) {
                if (e.cap > 0 && !seen[e.to]) {
                    seen[e.to]=true;
                    q.push(e.to);
                }
            }
        }
        return seen;
    }

    const DinicEdge& edge(int from, std::size_t index) const {
        return graph_[from][index];
    }

private:
    bool bfs(int s, int t) {
        std::fill(level_.begin(),level_.end(),-1);
        std::queue<int> q;
        level_[s]=0;
        q.push(s);
        while(!q.empty()) {
            const int v=q.front();
            q.pop();
            for (const auto& e : graph_[v]) {
                if (e.cap > 0 && level_[e.to] < 0) {
                    level_[e.to]=level_[v]+1;
                    q.push(e.to);
                }
            }
        }
        return level_[t] >= 0;
    }

    BigInt dfs(int v, int t, const BigInt& pushed) {
        if (pushed == 0 || v == t) return pushed;
        for (std::size_t& i=it_[v]; i<graph_[v].size(); ++i) {
            DinicEdge& e=graph_[v][i];
            if (e.cap <= 0 || level_[e.to] != level_[v]+1) continue;
            const BigInt next = e.cap < pushed ? e.cap : pushed;
            BigInt tr=dfs(e.to,t,next);
            if (tr == 0) continue;
            e.cap -= tr;
            graph_[e.to][e.rev].cap += tr;
            return tr;
        }
        return BigInt{0};
    }

    std::vector<std::vector<DinicEdge>> graph_;
    std::vector<int> level_;
    std::vector<std::size_t> it_;
};

void validate_input(const SeifertTurnFlowInput& input) {
    const std::size_t n=input.vertices.size();
    if (n == 0) throw std::invalid_argument("Seifert graph must contain at least one vertex");
    if (n > static_cast<std::size_t>(INT_MAX-4))
        throw std::invalid_argument("Seifert graph too large for current index backend");

    for (const auto& v : input.vertices) {
        if (v.positive_pieces < 0 || v.negative_pieces < 0)
            throw std::invalid_argument("Seifert piece counts must be nonnegative");
        if (v.positive_pieces + v.negative_pieces < 1)
            throw std::invalid_argument("every Seifert circle must contain at least one source piece");
        if (v.orientation != +1 && v.orientation != -1)
            throw std::invalid_argument("Seifert orientation must be +1 or -1");
    }

    for (const auto& e : input.edges) {
        if (e.u >= n || e.v >= n)
            throw std::invalid_argument("Seifert edge endpoint out of range");
    }
}

std::pair<std::size_t,std::size_t> canonical_endpoints(const SeifertEdge& e) {
    return e.u <= e.v ? std::pair{e.u,e.v} : std::pair{e.v,e.u};
}

std::size_t nonloop_edge_count(const SeifertTurnFlowInput& input) {
    std::size_t m=0;
    for (const auto& e : input.edges) if (e.u != e.v) ++m;
    return m;
}

struct BoundedArcRef {
    int from{};
    std::size_t edge_index{};
    BigInt lower{0};
    BigInt residual_initial{0};
    std::size_t input_edge_index{};
};

struct CirculationWitness {
    bool feasible{false};
    BigInt denominator{0};
    std::vector<BigInt> edge_x_scaled;
    std::vector<BigInt> source_t_scaled;
};

CirculationWitness solve_strict_circulation(const SeifertTurnFlowInput& input) {
    const std::size_t n=input.vertices.size();
    const std::size_t m=nonloop_edge_count(input);
    const std::size_t augmented_arc_count=m+n;
    const BigInt D = BigInt{2} * BigInt{augmented_arc_count+1};

    const int root=static_cast<int>(n);
    const int ss=static_cast<int>(n+1);
    const int tt=static_cast<int>(n+2);
    Dinic net(n+3);
    std::vector<BigInt> balance(n+1,BigInt{0});
    std::vector<BoundedArcRef> refs;
    refs.reserve(m);

    auto add_bounded = [&](int u, int v, const BigInt& lower, const BigInt& upper,
                           bool track, std::size_t input_index) {
        if (upper < lower) throw std::runtime_error("invalid shrunken circulation bounds");
        const BigInt residual=upper-lower;
        const std::size_t idx=net.add_edge(u,v,residual);
        balance[u]-=lower;
        balance[v]+=lower;
        if (track) refs.push_back(BoundedArcRef{u,idx,lower,residual,input_index});
    };

    std::vector<BigInt> q(n,BigInt{0});
    for (std::size_t i=0;i<input.edges.size();++i) {
        const auto [a,b]=canonical_endpoints(input.edges[i]);
        if (a == b) continue;
        q[a]+=1;
        q[b]-=1;
        add_bounded(static_cast<int>(a),static_cast<int>(b),
                    BigInt{1},BigInt{2}*D-BigInt{1},true,i);
    }

    for (std::size_t v=0;v<n;++v) {
        const auto& cap=input.vertices[v];
        const BigInt eps=cap.orientation;
        const BigInt L=BigInt{2}*eps-BigInt{2}*cap.positive_pieces+q[v];
        const BigInt U=BigInt{2}*eps+BigInt{2}*cap.negative_pieces+q[v];
        add_bounded(root,static_cast<int>(v),D*L+1,D*U-1,false,0);
    }

    BigInt total=0;
    for (std::size_t v=0;v<n+1;++v) {
        if (balance[v] > 0) {
            net.add_edge(ss,static_cast<int>(v),balance[v]);
            total += balance[v];
        } else if (balance[v] < 0) {
            net.add_edge(static_cast<int>(v),tt,-balance[v]);
        }
    }

    const BigInt pushed=net.max_flow(ss,tt,total);
    if (pushed != total) return {};

    CirculationWitness out;
    out.feasible=true;
    out.denominator=D;
    out.edge_x_scaled.assign(input.edges.size(),BigInt{0});
    out.source_t_scaled.assign(n,BigInt{0});

    for (const auto& ref : refs) {
        const BigInt residual_left=net.edge(ref.from,ref.edge_index).cap;
        const BigInt used=ref.residual_initial-residual_left;
        const BigInt g_scaled=ref.lower+used;
        out.edge_x_scaled[ref.input_edge_index]=g_scaled-D;
    }

    std::vector<BigInt> bx(n,BigInt{0});
    for (std::size_t i=0;i<input.edges.size();++i) {
        const auto [a,b]=canonical_endpoints(input.edges[i]);
        if (a == b) continue;
        const BigInt& x=out.edge_x_scaled[i];
        bx[a]+=x;
        bx[b]-=x;
    }

    for (std::size_t v=0;v<n;++v)
        out.source_t_scaled[v]=BigInt{2}*D*input.vertices[v].orientation-bx[v];

    return out;
}

BigInt abs_big(const BigInt& x) {
    return x < 0 ? -x : x;
}

struct CutCertificate {
    bool found{false};
    SeifertFlowProof proof{SeifertFlowProof::None};
    std::vector<bool> set;
};

bool exact_cut_violates(const SeifertTurnFlowInput& input,
                        const std::vector<bool>& set,
                        SeifertFlowProof proof) {
    if (set.size() != input.vertices.size()) return false;
    bool nonempty=false;
    BigInt P=0,N=0,eps=0,d=0;
    for (std::size_t v=0;v<set.size();++v) {
        if (!set[v]) continue;
        nonempty=true;
        P+=input.vertices[v].positive_pieces;
        N+=input.vertices[v].negative_pieces;
        eps+=input.vertices[v].orientation;
    }
    if (!nonempty) return false;

    for (const auto& e : input.edges) {
        if (e.u == e.v) continue;
        if (set[e.u] != set[e.v]) d+=1;
    }

    if (proof == SeifertFlowProof::PositiveCutObstruction)
        return BigInt{2}*eps >= BigInt{2}*P+d;
    if (proof == SeifertFlowProof::NegativeCutObstruction)
        return BigInt{2}*eps <= -BigInt{2}*N-d;
    return false;
}

CutCertificate find_weighted_cut_obstruction(const SeifertTurnFlowInput& input,
                                              SeifertFlowProof proof) {
    const std::size_t n=input.vertices.size();
    std::vector<BigInt> weight(n,BigInt{0});

    for (std::size_t v=0;v<n;++v) {
        if (proof == SeifertFlowProof::PositiveCutObstruction)
            weight[v]=BigInt{2}*(input.vertices[v].positive_pieces-input.vertices[v].orientation);
        else
            weight[v]=BigInt{2}*(input.vertices[v].negative_pieces+input.vertices[v].orientation);
    }

    BigInt finite_total=nonloop_edge_count(input);
    for (const auto& w : weight) finite_total+=abs_big(w);
    const BigInt INF=finite_total+1;

    for (std::size_t anchor=0;anchor<n;++anchor) {
        const int s=static_cast<int>(n);
        const int t=static_cast<int>(n+1);
        Dinic net(n+2);
        BigInt constant=0;
        BigInt limit=0;

        for (const auto& e : input.edges) {
            if (e.u == e.v) continue;
            net.add_edge(static_cast<int>(e.u),static_cast<int>(e.v),BigInt{1});
            net.add_edge(static_cast<int>(e.v),static_cast<int>(e.u),BigInt{1});
            limit+=1;
        }

        for (std::size_t v=0;v<n;++v) {
            if (weight[v] >= 0) {
                net.add_edge(static_cast<int>(v),t,weight[v]);
                limit+=weight[v];
            } else {
                net.add_edge(s,static_cast<int>(v),-weight[v]);
                limit+=-weight[v];
                constant+=weight[v];
            }
        }

        net.add_edge(s,static_cast<int>(anchor),INF);
        limit+=INF;

        const BigInt cut=net.max_flow(s,t,limit);
        const BigInt objective=cut+constant;
        if (objective <= 0) {
            const auto reachable=net.reachable(s);
            std::vector<bool> set(n,false);
            for (std::size_t v=0;v<n;++v) set[v]=reachable[v];
            if (set[anchor] && exact_cut_violates(input,set,proof))
                return CutCertificate{true,proof,std::move(set)};
        }
    }
    return {};
}

CutCertificate find_obstruction(const SeifertTurnFlowInput& input) {
    auto p=find_weighted_cut_obstruction(input,SeifertFlowProof::PositiveCutObstruction);
    if (p.found) return p;
    return find_weighted_cut_obstruction(input,SeifertFlowProof::NegativeCutObstruction);
}

SeifertFlowCertificate base_certificate(const SeifertTurnFlowInput& input) {
    SeifertFlowCertificate cert(input);
    cert.canonical_input_digest=sha256_hex(canonicalize_seifert_turn_flow_input(input));
    return cert;
}

SeifertFlowResult result_from_certificate(const SeifertFlowCertificate& cert) {
    SeifertFlowResult out(cert);
    out.status=cert.status;
    out.proof=cert.proof;
    out.assurance=cert.status == SeifertFlowStatus::Indeterminate
        ? ArithmeticAssurance::None : ArithmeticAssurance::Exact;
    out.termination=cert.status == SeifertFlowStatus::Indeterminate
        ? TerminationReason::BackendFailure : TerminationReason::Completed;
    return out;
}

bool verify_feasible_witness(const SeifertFlowCertificate& cert) {
    const auto& input=cert.input;
    const std::size_t n=input.vertices.size();
    if (cert.proof != SeifertFlowProof::InteriorFlowWitness) return false;
    if (cert.denominator <= 0) return false;
    if (cert.scaled_edge_flow.size() != input.edges.size()) return false;
    if (cert.scaled_source_turn.size() != n) return false;
    if (!cert.cut_vertices.empty()) return false;

    const BigInt& D=cert.denominator;
    std::vector<BigInt> bx(n,BigInt{0});
    for (std::size_t i=0;i<input.edges.size();++i) {
        const auto [a,b]=canonical_endpoints(input.edges[i]);
        const BigInt& x=cert.scaled_edge_flow[i];
        if (a == b) {
            if (x != 0) return false;
            continue;
        }
        if (!(-D < x && x < D)) return false;
        bx[a]+=x;
        bx[b]-=x;
    }

    for (std::size_t v=0;v<n;++v) {
        const BigInt& t=cert.scaled_source_turn[v];
        const auto& cap=input.vertices[v];
        if (!(-BigInt{2}*D*cap.negative_pieces < t &&
              t < BigInt{2}*D*cap.positive_pieces))
            return false;
        if (t+bx[v] != BigInt{2}*D*cap.orientation)
            return false;
    }
    return true;
}

bool metadata_ok(const SeifertFlowCertificate& cert) {
    return cert.schema == kSchema &&
           cert.schema_version == kSchemaVersion &&
           cert.theorem_id == kTheoremId &&
           cert.source_digest.algorithm == "SHA-256" &&
           cert.source_digest.value == kSourceDigest &&
           cert.canonical_input_digest == sha256_hex(canonicalize_seifert_turn_flow_input(cert.input));
}

SeifertFlowResult invalid_certificate(const SeifertFlowCertificate& cert) {
    SeifertFlowResult out(cert);
    out.status=SeifertFlowStatus::Indeterminate;
    out.proof=SeifertFlowProof::None;
    out.assurance=ArithmeticAssurance::None;
    out.termination=TerminationReason::BackendFailure;
    return out;
}

} // namespace

std::string canonicalize_seifert_turn_flow_input(const SeifertTurnFlowInput& input) {
    validate_input(input);
    std::vector<std::pair<std::size_t,std::size_t>> edges;
    edges.reserve(input.edges.size());
    for (const auto& e : input.edges) edges.push_back(canonical_endpoints(e));
    std::sort(edges.begin(),edges.end());

    std::ostringstream os;
    os << "SEIFERT-TURN-FLOW|vertices=" << input.vertices.size() << "|caps=";
    for (std::size_t i=0;i<input.vertices.size();++i) {
        if (i) os << ";";
        const auto& v=input.vertices[i];
        os << i << ":" << v.positive_pieces << "/" << v.negative_pieces << "/" << v.orientation;
    }
    os << "|edges=";
    for (std::size_t i=0;i<edges.size();++i) {
        if (i) os << ";";
        os << edges[i].first << "-" << edges[i].second;
    }
    return os.str();
}

std::string serialize_seifert_turn_flow_certificate(const SeifertFlowCertificate& c) {
    std::ostringstream os;
    os << "{\"schema\":\"" << c.schema
       << "\",\"schema_version\":\"" << c.schema_version
       << "\",\"theorem_id\":\"" << c.theorem_id
       << "\",\"source_digest\":\"" << c.source_digest.value
       << "\",\"canonical_input_digest\":\"" << c.canonical_input_digest
       << "\",\"status\":" << static_cast<int>(c.status)
       << ",\"proof\":" << static_cast<int>(c.proof)
       << ",\"denominator\":\"" << c.denominator << "\""
       << ",\"edge_flow\":[";
    for (std::size_t i=0;i<c.scaled_edge_flow.size();++i) {
        if (i) os << ",";
        os << "\"" << c.scaled_edge_flow[i] << "\"";
    }
    os << "],\"source_turn\":[";
    for (std::size_t i=0;i<c.scaled_source_turn.size();++i) {
        if (i) os << ",";
        os << "\"" << c.scaled_source_turn[i] << "\"";
    }
    os << "],\"cut\":[";
    bool first=true;
    for (std::size_t i=0;i<c.cut_vertices.size();++i) {
        if (!c.cut_vertices[i]) continue;
        if (!first) os << ",";
        first=false;
        os << i;
    }
    os << "]}";
    return os.str();
}

SeifertFlowResult solve_seifert_turn_flow(const SeifertTurnFlowInput& input) {
    validate_input(input);
    auto cert=base_certificate(input);

    try {
        auto witness=solve_strict_circulation(input);
        if (witness.feasible) {
            cert.status=SeifertFlowStatus::Feasible;
            cert.proof=SeifertFlowProof::InteriorFlowWitness;
            cert.denominator=witness.denominator;
            cert.scaled_edge_flow=std::move(witness.edge_x_scaled);
            cert.scaled_source_turn=std::move(witness.source_t_scaled);
            auto verified=verify_seifert_turn_flow_certificate(cert);
            if (verified.status != SeifertFlowStatus::Feasible)
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        auto obstruction=find_obstruction(input);
        if (obstruction.found) {
            cert.status=SeifertFlowStatus::Infeasible;
            cert.proof=obstruction.proof;
            cert.cut_vertices=std::move(obstruction.set);
            auto verified=verify_seifert_turn_flow_certificate(cert);
            if (verified.status != SeifertFlowStatus::Infeasible)
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }
    } catch (...) {
        return invalid_certificate(cert);
    }
    return invalid_certificate(cert);
}

SeifertFlowResult verify_seifert_turn_flow_certificate(const SeifertFlowCertificate& cert) {
    try {
        validate_input(cert.input);
        if (!metadata_ok(cert)) return invalid_certificate(cert);

        if (cert.status == SeifertFlowStatus::Feasible) {
            if (!verify_feasible_witness(cert)) return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        if (cert.status == SeifertFlowStatus::Infeasible) {
            if (cert.denominator != 0 ||
                !cert.scaled_edge_flow.empty() ||
                !cert.scaled_source_turn.empty())
                return invalid_certificate(cert);
            if (!exact_cut_violates(cert.input,cert.cut_vertices,cert.proof))
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        if (cert.proof != SeifertFlowProof::None)
            return invalid_certificate(cert);
        return result_from_certificate(cert);
    } catch (...) {
        return invalid_certificate(cert);
    }
}

} // namespace pcg
