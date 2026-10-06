#include <pcg/post_closure.hpp>

#include <algorithm>
#include <array>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace pcg {
namespace {

constexpr std::size_t kNoIndex=std::numeric_limits<std::size_t>::max();
const char* kSchema="pcg-rigid-post-closure-certificate";
const char* kSchemaVersion="1.0";
const char* kTheoremId="PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001";
const char* kSourceDigest="625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208";

std::string closure_binding(const NetworkClosureCertificate& c) {
    return sha256_hex(
        canonicalize_network_closure_input(c.input)+"|"+
        serialize_network_closure_certificate(c));
}

RigidPostClosureCertificate base_certificate(const NetworkClosureCertificate& c) {
    RigidPostClosureCertificate out(c);
    out.closure_binding_digest=closure_binding(c);
    return out;
}

RigidPostClosureResult result_from_certificate(const RigidPostClosureCertificate& c) {
    RigidPostClosureResult out(c);
    out.decision=c.decision;
    out.proof=c.proof;
    out.assurance=c.assurance;
    out.termination=c.termination;
    return out;
}

RigidPostClosureResult invalid_certificate(const RigidPostClosureCertificate& c) {
    RigidPostClosureCertificate bad=c;
    bad.decision=Decision::Indeterminate;
    bad.proof=RigidPostClosureProof::BackendFailure;
    bad.assurance=ArithmeticAssurance::None;
    bad.termination=TerminationReason::BackendFailure;
    return result_from_certificate(bad);
}

bool metadata_ok(const RigidPostClosureCertificate& c) {
    return c.schema==kSchema &&
           c.schema_version==kSchemaVersion &&
           c.theorem_id==kTheoremId &&
           c.source_digest.algorithm=="SHA-256" &&
           c.source_digest.value==kSourceDigest &&
           c.closure_binding_digest==closure_binding(c.closure_certificate);
}

bool algebraic_equal(const AlgebraicChordValue& a,const AlgebraicChordValue& b) {
    return a.coefficients==b.coefficients;
}

AlgebraicChordValue algebraic_subtract(const AlgebraicChordValue& a,
                                       const AlgebraicChordValue& b) {
    if (a.coefficients.size()!=b.coefficients.size())
        throw std::invalid_argument("algebraic basis mismatch");
    AlgebraicChordValue out;
    out.coefficients.resize(a.coefficients.size(),Rational(BigInt{0}));
    for (std::size_t i=0;i<a.coefficients.size();++i)
        out.coefficients[i]=a.coefficients[i]-b.coefficients[i];
    return out;
}

enum class VertexRelation { Distinct, Coincident, Indeterminate };

VertexRelation vertex_relation(unsigned order,
                               const AlgebraicPoint2& a,
                               const AlgebraicPoint2& b) {
    if (algebraic_equal(a.x,b.x) && algebraic_equal(a.y,b.y))
        return VertexRelation::Coincident;

    const auto dx=certified_algebraic_real(order,algebraic_subtract(a.x,b.x));
    if (dx.lower()>0.0 || dx.upper()<0.0) return VertexRelation::Distinct;
    const auto dy=certified_algebraic_real(order,algebraic_subtract(a.y,b.y));
    if (dy.lower()>0.0 || dy.upper()<0.0) return VertexRelation::Distinct;
    return VertexRelation::Indeterminate;
}

CertifiedPoint2 certified_point(unsigned order,const AlgebraicPoint2& p) {
    return {certified_algebraic_real(order,p.x),
            certified_algebraic_real(order,p.y)};
}

std::vector<std::size_t> shared_vertices(const NetworkClosureInput& in,
                                         std::size_t a,std::size_t b) {
    const std::size_t m=in.turns.size();
    const std::array<std::size_t,2> ea{
        in.trace_vertices[a],in.trace_vertices[(a+1)%m]};
    const std::array<std::size_t,2> eb{
        in.trace_vertices[b],in.trace_vertices[(b+1)%m]};
    std::vector<std::size_t> out;
    for (auto x:ea) for (auto y:eb) {
        if (x!=y) continue;
        if (std::find(out.begin(),out.end(),x)==out.end()) out.push_back(x);
    }
    return out;
}

SharedStartArc shared_at_vertex(const NetworkClosureInput& in,
                                const ProjectivelyRigidEmbedding& emb,
                                const std::vector<PositiveInterval>& chords,
                                std::size_t e,std::size_t v) {
    const std::size_t m=in.turns.size();
    const std::size_t tail=in.trace_vertices[e];
    const std::size_t head=in.trace_vertices[(e+1)%m];
    if (tail==v)
        return SharedStartArc{emb.tangent_phase_pi[e],in.turns[e],chords[e]};
    if (head==v)
        return SharedStartArc{
            emb.tangent_phase_pi[e]+in.turns[e].pi+PiRational{1,1},
            Turn{-in.turns[e].pi},
            chords[e]};
    throw std::invalid_argument("edge is not incident to requested vertex");
}

RigidPostClosureCertificate evaluate(const NetworkClosureCertificate& closure) {
    auto cert=base_certificate(closure);

    const auto verified=verify_network_closure_certificate(closure);
    if (verified.termination==TerminationReason::BackendFailure) {
        cert.proof=RigidPostClosureProof::BackendFailure;
        return cert;
    }
    if (verified.status==NetworkClosureStatus::Indeterminate) {
        cert.decision=Decision::Indeterminate;
        cert.proof=RigidPostClosureProof::PairCertificationLimit;
        cert.termination=verified.termination;
        return cert;
    }
    if (verified.status==NetworkClosureStatus::NotClosed) {
        cert.decision=Decision::NotRealizable;
        cert.proof=RigidPostClosureProof::NetworkClosureObstruction;
        cert.assurance=verified.assurance;
        cert.termination=TerminationReason::Completed;
        return cert;
    }

    const auto reconstructed=reconstruct_projectively_rigid_embedding(closure);
    if (reconstructed.status==ProjectivelyRigidEmbeddingStatus::NotProjectivelyRigid) {
        cert.decision=Decision::Indeterminate;
        cert.proof=RigidPostClosureProof::HigherDimensionalClosure;
        cert.termination=TerminationReason::Completed;
        return cert;
    }
    if (reconstructed.status!=ProjectivelyRigidEmbeddingStatus::Ready) {
        cert.decision=Decision::Indeterminate;
        cert.proof=reconstructed.status==ProjectivelyRigidEmbeddingStatus::InvalidCertificate
            ? RigidPostClosureProof::BackendFailure
            : RigidPostClosureProof::PairCertificationLimit;
        cert.termination=reconstructed.status==ProjectivelyRigidEmbeddingStatus::InvalidCertificate
            ? TerminationReason::BackendFailure
            : TerminationReason::PrecisionLimit;
        return cert;
    }

    const auto& in=closure.input;
    const auto& emb=reconstructed.embedding;
    const std::size_t m=in.turns.size();

    // A Gauss double point may occur twice. Higher multiplicity belongs to a
    // different singularity model and is outside this Phase 6A specialization.
    std::vector<std::vector<std::size_t>> occurrences(in.vertex_count);
    for (std::size_t e=0;e<m;++e)
        occurrences[in.trace_vertices[e]].push_back(e);
    for (std::size_t v=0;v<occurrences.size();++v) {
        if (occurrences[v].size()>2) {
            cert.decision=Decision::Unsupported;
            cert.proof=RigidPostClosureProof::UnsupportedTraceMultiplicity;
            cert.assurance=ArithmeticAssurance::Exact;
            cert.termination=TerminationReason::Completed;
            cert.witness_a=v;
            return cert;
        }
        if (occurrences[v].size()==2) {
            const auto delta=emb.tangent_phase_pi[occurrences[v][1]]
                           - emb.tangent_phase_pi[occurrences[v][0]];
            if (sin_pi_sign(delta)==0) {
                cert.decision=Decision::Unsupported;
                cert.proof=RigidPostClosureProof::NonTransversePrescribedPassage;
                cert.assurance=ArithmeticAssurance::Exact;
                cert.termination=TerminationReason::Completed;
                cert.witness_a=v;
                return cert;
            }
        }
    }

    for (std::size_t v=0;v<emb.vertices.size();++v) {
        for (std::size_t w=v+1;w<emb.vertices.size();++w) {
            const auto rel=vertex_relation(emb.cyclotomic_order,
                                           emb.vertices[v],emb.vertices[w]);
            if (rel==VertexRelation::Coincident) {
                cert.decision=Decision::NotRealizable;
                cert.proof=RigidPostClosureProof::ExactVertexCollision;
                cert.assurance=ArithmeticAssurance::Exact;
                cert.termination=TerminationReason::Completed;
                cert.witness_a=v;
                cert.witness_b=w;
                return cert;
            }
            if (rel==VertexRelation::Indeterminate) {
                cert.decision=Decision::Indeterminate;
                cert.proof=RigidPostClosureProof::PairCertificationLimit;
                cert.termination=TerminationReason::PrecisionLimit;
                cert.witness_a=v;
                cert.witness_b=w;
                return cert;
            }
        }
    }

    std::vector<CertifiedPoint2> vertices;
    vertices.reserve(emb.vertices.size());
    for (const auto& p:emb.vertices)
        vertices.push_back(certified_point(emb.cyclotomic_order,p));

    std::vector<PositiveInterval> chords;
    chords.reserve(m);
    for (const auto& c:emb.chord_magnitudes) {
        const auto ci=certified_algebraic_real(emb.cyclotomic_order,c);
        if (!(ci.lower()>0.0)) {
            cert.decision=Decision::Indeterminate;
            cert.proof=RigidPostClosureProof::PairCertificationLimit;
            cert.termination=TerminationReason::PrecisionLimit;
            return cert;
        }
        chords.emplace_back(ci);
    }

    std::vector<CertifiedFiniteArc> arcs;
    arcs.reserve(m);
    for (std::size_t e=0;e<m;++e) {
        const std::size_t tail=in.trace_vertices[e];
        const std::size_t head=in.trace_vertices[(e+1)%m];
        arcs.push_back(CertifiedFiniteArc{
            vertices[tail],vertices[head],emb.tangent_phase_pi[e],
            in.turns[e],chords[e]});
    }

    for (std::size_t e=0;e<m;++e) {
        for (std::size_t f=e+1;f<m;++f) {
            ++cert.checked_pairs;
            const auto shared=shared_vertices(in,e,f);
            CertifiedTruth hit=CertifiedTruth::Indeterminate;

            if (shared.empty()) {
                hit=pcg_hit(arcs[e],arcs[f]);
            } else if (shared.size()==1) {
                const auto a=shared_at_vertex(in,emb,chords,e,shared.front());
                const auto b=shared_at_vertex(in,emb,chords,f,shared.front());

                // If the two outgoing tangent lines agree exactly at the
                // prescribed common endpoint, distinct supporting circles are
                // tangent there. Two distinct circles cannot then have a
                // second common point, so the pair is certified clear without
                // asking the generic ExtraHit interval construction to
                // separate the tangent boundary.
                const PiRational tangent_delta=
                    b.tangent_phase_pi-a.tangent_phase_pi;
                if (sin_pi_sign(tangent_delta)==0) {
                    const auto relation=shared_start_support_relation(a,b);
                    if (relation==SupportRelation::Distinct)
                        hit=CertifiedTruth::False;
                    else
                        hit=CertifiedTruth::Indeterminate;
                } else {
                    hit=pcg_extra_hit(a,b);
                }
            } else if (shared.size()==2) {
                const auto a=shared_at_vertex(in,emb,chords,e,shared.front());
                const auto b=shared_at_vertex(in,emb,chords,f,shared.front());
                const auto relation=shared_start_support_relation(a,b);
                if (relation==SupportRelation::Distinct)
                    hit=CertifiedTruth::False;
                else
                    hit=CertifiedTruth::Indeterminate;
            } else {
                cert.decision=Decision::Unsupported;
                cert.proof=RigidPostClosureProof::UnsupportedTraceMultiplicity;
                cert.assurance=ArithmeticAssurance::Exact;
                cert.termination=TerminationReason::Completed;
                cert.witness_a=e;
                cert.witness_b=f;
                return cert;
            }

            if (hit==CertifiedTruth::True) {
                cert.decision=Decision::NotRealizable;
                cert.proof=RigidPostClosureProof::CertifiedUnintendedIntersection;
                cert.assurance=ArithmeticAssurance::CertifiedNumerical;
                cert.termination=TerminationReason::Completed;
                cert.witness_a=e;
                cert.witness_b=f;
                return cert;
            }
            if (hit==CertifiedTruth::Indeterminate) {
                cert.decision=Decision::Indeterminate;
                cert.proof=RigidPostClosureProof::PairCertificationLimit;
                cert.termination=TerminationReason::PrecisionLimit;
                cert.witness_a=e;
                cert.witness_b=f;
                return cert;
            }
        }
    }

    cert.decision=Decision::Realizable;
    cert.proof=RigidPostClosureProof::AllPairsClear;
    cert.assurance=ArithmeticAssurance::CertifiedNumerical;
    cert.termination=TerminationReason::Completed;
    cert.witness_a=kNoIndex;
    cert.witness_b=kNoIndex;
    return cert;
}

} // namespace

std::string serialize_rigid_post_closure_certificate(
    const RigidPostClosureCertificate& c) {
    std::ostringstream os;
    os << "{\"schema\":\"" << c.schema
       << "\",\"schema_version\":\"" << c.schema_version
       << "\",\"theorem_id\":\"" << c.theorem_id
       << "\",\"source_digest\":\"" << c.source_digest.value
       << "\",\"closure_binding_digest\":\"" << c.closure_binding_digest
       << "\",\"decision\":" << static_cast<int>(c.decision)
       << ",\"proof\":" << static_cast<int>(c.proof)
       << ",\"assurance\":" << static_cast<int>(c.assurance)
       << ",\"termination\":" << static_cast<int>(c.termination)
       << ",\"witness_a\":" << c.witness_a
       << ",\"witness_b\":" << c.witness_b
       << ",\"checked_pairs\":" << c.checked_pairs
       << "}";
    return os.str();
}

RigidPostClosureResult solve_projectively_rigid_post_closure(
    const NetworkClosureCertificate& closure_certificate) {
    return result_from_certificate(evaluate(closure_certificate));
}

RigidPostClosureResult verify_projectively_rigid_post_closure_certificate(
    const RigidPostClosureCertificate& cert) {
    try {
        if (!metadata_ok(cert)) return invalid_certificate(cert);
        const auto expected=evaluate(cert.closure_certificate);
        if (expected.decision!=cert.decision ||
            expected.proof!=cert.proof ||
            expected.assurance!=cert.assurance ||
            expected.termination!=cert.termination ||
            expected.witness_a!=cert.witness_a ||
            expected.witness_b!=cert.witness_b ||
            expected.checked_pairs!=cert.checked_pairs)
            return invalid_certificate(cert);
        return result_from_certificate(cert);
    } catch (...) {
        return invalid_certificate(cert);
    }
}

} // namespace pcg
