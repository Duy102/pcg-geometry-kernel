#include <pcg/abcabc.hpp>
#include <array>
#include <utility>
#include <sstream>

namespace pcg {
namespace {

constexpr ABCABCSignPattern kR0Representative{{+1,+1,+1,-1,-1,-1}};
constexpr std::array<ABCABCSignPattern,4> kR2Representatives{{
    ABCABCSignPattern{{+1,+1,+1,+1,+1,+1}},
    ABCABCSignPattern{{+1,+1,+1,+1,+1,-1}},
    ABCABCSignPattern{{+1,+1,+1,-1,+1,-1}},
    ABCABCSignPattern{{+1,-1,+1,-1,+1,-1}}
}};

ABCABCSignPattern rotate_left(const ABCABCSignPattern& s, std::size_t shift) {
    ABCABCSignPattern out{};
    for (std::size_t i=0;i<6;++i) out[i]=s[(i+shift)%6];
    return out;
}

ABCABCSignPattern reversed_signs(const ABCABCSignPattern& s) {
    ABCABCSignPattern out{};
    for (std::size_t i=0;i<6;++i) out[i]=-s[i];
    return out;
}

bool is_valid_sign_pattern(const ABCABCSignPattern& s) {
    for (int x : s) if (x != +1 && x != -1) return false;
    return true;
}

bool in_cyclic_orbit(const ABCABCSignPattern& s, const ABCABCSignPattern& representative) {
    for (std::size_t shift=0;shift<6;++shift)
        if (s == rotate_left(representative,shift)) return true;
    return false;
}

bool r2_class_realizable(const ABCABCSignPattern& s) {
    for (const auto& representative : kR2Representatives)
        if (in_cyclic_orbit(s,representative)) return true;
    return false;
}

PiRational half(PiRational q) { return q / 2; }

std::array<PiRational,6> qturns(const ABCABCInput& in) {
    std::array<PiRational,6> q{};
    for (std::size_t i=0;i<6;++i) q[i]=in.turns[i].pi;
    return q;
}

bool phase_ok(const std::array<PiRational,6>& q) {
    const PiRational p1 = half(q[0]) + q[1] + q[2] + half(q[3]);
    const PiRational p2 = half(q[1]) + q[2] + q[3] + half(q[4]);
    const PiRational p3 = half(q[2]) + q[3] + q[4] + half(q[5]);
    return is_even_integer(p1) && is_even_integer(p2) && is_even_integer(p3);
}

std::array<PiRational,6> phases(const std::array<PiRational,6>& q) {
    std::array<PiRational,6> p{};
    p[0]=PiRational{0,1};
    for(std::size_t i=1;i<6;++i) p[i]=p[i-1]+q[i-1];
    return p;
}

std::array<PiRational,3> alphas(const std::array<PiRational,6>& q, const std::array<PiRational,6>& p) {
    return {p[0]+half(q[0]), p[1]+half(q[1]), p[2]+half(q[2])};
}

bool positive_closure_ok(const std::array<PiRational,3>& a) {
    const int s1=sin_pi_sign(a[2]-a[1]);
    const int s2=sin_pi_sign(a[0]-a[2]);
    const int s3=sin_pi_sign(a[1]-a[0]);
    return s1!=0 && s1==s2 && s2==s3;
}

PositiveInterval chord_from_sine(PiRational q) {
    const auto x=certified_abs(certified_sin_pi(q));
    return PositiveInterval(x);
}

SharedStartArc forward_arc(std::size_t i,
                           const std::array<PiRational,6>& p,
                           const ABCABCInput& in,
                           const std::array<PositiveInterval,6>& chords) {
    return SharedStartArc{p[i], in.turns[i], chords[i]};
}

SharedStartArc reverse_arc(std::size_t i,
                           const std::array<PiRational,6>& p,
                           const ABCABCInput& in,
                           const std::array<PositiveInterval,6>& chords) {
    return SharedStartArc{p[i]+in.turns[i].pi+PiRational{1,1}, Turn{-in.turns[i].pi}, chords[i]};
}

using ArcPair = std::pair<SharedStartArc,SharedStartArc>;

std::array<ArcPair,3> paired_endpoint_pairs(const std::array<PiRational,6>& p,
                                            const ABCABCInput& in,
                                            const std::array<PositiveInterval,6>& c) {
    return {{
        {forward_arc(0,p,in,c), forward_arc(3,p,in,c)},
        {forward_arc(1,p,in,c), forward_arc(4,p,in,c)},
        {forward_arc(2,p,in,c), forward_arc(5,p,in,c)}
    }};
}

std::array<ArcPair,6> source_adjacent_pairs(const std::array<PiRational,6>& p,
                                            const ABCABCInput& in,
                                            const std::array<PositiveInterval,6>& c) {
    return {{
        {reverse_arc(0,p,in,c), forward_arc(1,p,in,c)},
        {reverse_arc(1,p,in,c), forward_arc(2,p,in,c)},
        {reverse_arc(2,p,in,c), forward_arc(3,p,in,c)},
        {reverse_arc(3,p,in,c), forward_arc(4,p,in,c)},
        {reverse_arc(4,p,in,c), forward_arc(5,p,in,c)},
        {reverse_arc(5,p,in,c), forward_arc(0,p,in,c)}
    }};
}

std::array<ArcPair,6> cross_passage_pairs(const std::array<PiRational,6>& p,
                                          const ABCABCInput& in,
                                          const std::array<PositiveInterval,6>& c) {
    return {{
        {forward_arc(0,p,in,c), reverse_arc(2,p,in,c)},
        {forward_arc(3,p,in,c), reverse_arc(5,p,in,c)},
        {forward_arc(1,p,in,c), reverse_arc(3,p,in,c)},
        {forward_arc(4,p,in,c), reverse_arc(0,p,in,c)},
        {forward_arc(2,p,in,c), reverse_arc(4,p,in,c)},
        {forward_arc(5,p,in,c), reverse_arc(1,p,in,c)}
    }};
}

CertifiedTruth non_antipodal_second_for_arc(PiRational relative_phase,
                                             const SharedStartIntersection& inter) {
    const auto s = certified_sin_pi(relative_phase).interval();
    const auto c = certified_cos_pi(relative_phase).interval();
    const auto& qx = inter.qx.interval();
    const auto& qy = inter.qy.interval();
    const auto xr = qx*c + qy*s;
    const auto yr = -qx*s + qy*c;

    // Antipodal means the second intersection is exactly on the negative
    // initial-tangent ray. Prove non-antipodal by separating either x or y.
    if (xr.lower() >= 0.0) return CertifiedTruth::True;
    if (yr.lower() > 0.0 || yr.upper() < 0.0) return CertifiedTruth::True;
    if (xr.upper() < 0.0 && yr.lower() == 0.0 && yr.upper() == 0.0)
        return CertifiedTruth::False;
    return CertifiedTruth::Indeterminate;
}

ABCABCGenericityStatus genericity_from_precomputed(const std::array<PiRational,6>& p,
                                                    const ABCABCInput& in,
                                                    const std::array<PositiveInterval,6>& c) {
    const auto paired = paired_endpoint_pairs(p,in,c);
    const auto adjacent = source_adjacent_pairs(p,in,c);
    const auto cross = cross_passage_pairs(p,in,c);

    bool unresolved=false;
    auto inspect_support = [&](const auto& pairs) {
        for (const auto& [x,y] : pairs) {
            const auto relation = shared_start_support_relation(x,y);
            if (relation == SupportRelation::Coincident)
                return ABCABCGenericityStatus::Violated;
            if (relation == SupportRelation::Indeterminate ||
                relation == SupportRelation::BackendFailure)
                unresolved=true;
        }
        return ABCABCGenericityStatus::Satisfied;
    };

    if (inspect_support(paired) == ABCABCGenericityStatus::Violated ||
        inspect_support(adjacent) == ABCABCGenericityStatus::Violated ||
        inspect_support(cross) == ABCABCGenericityStatus::Violated)
        return ABCABCGenericityStatus::Violated;

    // Cross-passage arcs meet at prescribed double points, so their common
    // start must be transverse (not tangent), and the second intersection
    // must avoid the principal-argument antipodal branch cut.
    for (const auto& [x,y] : cross) {
        const auto inter = shared_start_second_intersection(x,y);
        if (inter.status == GeometricStatus::TangentAtStart ||
            inter.status == GeometricStatus::CoincidentSupport)
            return ABCABCGenericityStatus::Violated;
        if (inter.status != GeometricStatus::RegularSecondIntersection) {
            unresolved=true;
            continue;
        }
        const PiRational delta = y.tangent_phase_pi - x.tangent_phase_pi;
        const auto na = non_antipodal_second_for_arc(PiRational{0,1},inter);
        const auto nb = non_antipodal_second_for_arc(delta,inter);
        if (na == CertifiedTruth::False || nb == CertifiedTruth::False)
            return ABCABCGenericityStatus::Violated;
        if (na == CertifiedTruth::Indeterminate || nb == CertifiedTruth::Indeterminate)
            unresolved=true;
    }

    return unresolved ? ABCABCGenericityStatus::Indeterminate
                      : ABCABCGenericityStatus::Satisfied;
}

struct Evaluation {
    Decision decision;
    ProofKind kind;
    ArithmeticAssurance assurance;
    ABCABCProofReason reason;
};

Evaluation evaluate(const ABCABCInput& in) {
    const auto q=qturns(in);
    if (!phase_ok(q)) return {Decision::NotRealizable, ProofKind::TheoremInstantiation, ArithmeticAssurance::Exact, ABCABCProofReason::PhaseCongruenceObstruction};
    const auto p=phases(q);
    const auto a=alphas(q,p);
    if (!positive_closure_ok(a)) return {Decision::NotRealizable, ProofKind::TheoremInstantiation, ArithmeticAssurance::Exact, ABCABCProofReason::PositiveClosureObstruction};

    std::array<PositiveInterval,6> c{
        chord_from_sine(a[2]-a[1]),
        chord_from_sine(a[0]-a[2]),
        chord_from_sine(a[1]-a[0]),
        chord_from_sine(a[2]-a[1]),
        chord_from_sine(a[0]-a[2]),
        chord_from_sine(a[1]-a[0])
    };

    const auto genericity = genericity_from_precomputed(p,in,c);
    if (genericity == ABCABCGenericityStatus::Violated)
        return {Decision::Unsupported, ProofKind::FinitePredicateCertificate,
                ArithmeticAssurance::CertifiedNumerical, ABCABCProofReason::TheoremDomainViolation};
    if (genericity == ABCABCGenericityStatus::Indeterminate)
        return {Decision::Indeterminate, ProofKind::None,
                ArithmeticAssurance::None, ABCABCProofReason::NumericalIndeterminacy};

    const auto pairs = cross_passage_pairs(p,in,c);

    bool indeterminate=false;
    for (const auto& [x,y] : pairs) {
        const auto hit=pcg_extra_hit(x,y);
        if (hit==CertifiedTruth::True)
            return {Decision::NotRealizable, ProofKind::FinitePredicateCertificate, ArithmeticAssurance::CertifiedNumerical, ABCABCProofReason::ExtraHitObstruction};
        if (hit==CertifiedTruth::Indeterminate) indeterminate=true;
    }
    if (indeterminate)
        return {Decision::Indeterminate, ProofKind::None, ArithmeticAssurance::None, ABCABCProofReason::NumericalIndeterminacy};
    return {Decision::Realizable, ProofKind::CompositeCertificate, ArithmeticAssurance::CertifiedNumerical, ABCABCProofReason::AllConditionsSatisfied};
}

ABCABCResult result_from_eval(const ABCABCInput& in, const Evaluation& e) {
    ABCABCCertificate cert(in);
    cert.canonical_input_digest=sha256_hex(canonicalize_abcabc_input(in));
    cert.reason=e.reason;
    cert.proof_kind=e.kind;
    cert.assurance=e.assurance;
    return ABCABCResult{e.decision,e.kind,e.assurance,
        e.decision==Decision::Indeterminate ? TerminationReason::PrecisionLimit : TerminationReason::Completed,
        cert};
}
}

ABCABCGenericityStatus certify_abcabc_genericity(const ABCABCInput& input) {
    try {
        const auto q=qturns(input);
        if (!phase_ok(q)) return ABCABCGenericityStatus::NotApplicable;
        const auto p=phases(q);
        const auto a=alphas(q,p);
        if (!positive_closure_ok(a)) return ABCABCGenericityStatus::NotApplicable;
        std::array<PositiveInterval,6> c{
            chord_from_sine(a[2]-a[1]),
            chord_from_sine(a[0]-a[2]),
            chord_from_sine(a[1]-a[0]),
            chord_from_sine(a[2]-a[1]),
            chord_from_sine(a[0]-a[2]),
            chord_from_sine(a[1]-a[0])
        };
        return genericity_from_precomputed(p,input,c);
    } catch (...) {
        return ABCABCGenericityStatus::Indeterminate;
    }
}

Decision classify_abcabc_sign_rotation(const ABCABCSignPattern& signs, int rotation) {
    if (!is_valid_sign_pattern(signs))
        throw std::invalid_argument("ABCABC sign pattern entries must be +1 or -1");

    if (rotation == 0)
        return in_cyclic_orbit(signs,kR0Representative)
            ? Decision::Realizable : Decision::NotRealizable;

    if (rotation == 2)
        return r2_class_realizable(signs)
            ? Decision::Realizable : Decision::NotRealizable;

    if (rotation == -2)
        return r2_class_realizable(reversed_signs(signs))
            ? Decision::Realizable : Decision::NotRealizable;

    return Decision::NotRealizable;
}

std::string canonicalize_abcabc_input(const ABCABCInput& input) {
    std::ostringstream os;
    os << "ABCABC|turns=";
    for (std::size_t i=0;i<input.turns.size();++i) {
        if (i) os << ",";
        os << rational_string(input.turns[i].pi);
    }
    return os.str();
}

std::string serialize_abcabc_certificate(const ABCABCCertificate& c) {
    std::ostringstream os;
    os << "{\"schema\":\"" << c.schema
       << "\",\"schema_version\":\"" << c.schema_version
       << "\",\"theorem_id\":\"" << c.theorem_id
       << "\",\"source_digest_algorithm\":\"" << c.source_digest.algorithm
       << "\",\"source_digest\":\"" << c.source_digest.value
       << "\",\"canonical_input_digest\":\"" << c.canonical_input_digest
       << "\",\"turns_pi\":[";
    for (std::size_t i=0;i<c.input.turns.size();++i) {
        if (i) os << ",";
        os << "\"" << rational_string(c.input.turns[i].pi) << "\"";
    }
    os << "],\"reason\":" << static_cast<int>(c.reason)
       << ",\"proof_kind\":" << static_cast<int>(c.proof_kind)
       << ",\"arithmetic_assurance\":" << static_cast<int>(c.assurance) << "}";
    return os.str();
}

ABCABCResult solve_abcabc(const ABCABCInput& input) {
    return result_from_eval(input,evaluate(input));
}

} // namespace pcg
