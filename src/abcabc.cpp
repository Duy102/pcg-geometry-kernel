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

    const std::array<std::pair<SharedStartArc,SharedStartArc>,6> pairs{{
        {forward_arc(0,p,in,c), reverse_arc(2,p,in,c)},
        {forward_arc(3,p,in,c), reverse_arc(5,p,in,c)},
        {forward_arc(1,p,in,c), reverse_arc(3,p,in,c)},
        {forward_arc(4,p,in,c), reverse_arc(0,p,in,c)},
        {forward_arc(2,p,in,c), reverse_arc(4,p,in,c)},
        {forward_arc(5,p,in,c), reverse_arc(1,p,in,c)}
    }};

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

ABCABCResult verify_abcabc_certificate(const ABCABCCertificate& certificate) {
    if (certificate.schema != "pcg-abcabc-certificate" || certificate.schema_version != "1.0" ||
        certificate.theorem_id != "PCG-ABCABC-THM-001" ||
        certificate.source_digest.algorithm != "SHA-256" ||
        certificate.source_digest.value != "2aa832a00031035d964c0d9362e4409861e19b97a9d764f3123a42edef400beb" ||
        certificate.canonical_input_digest != sha256_hex(canonicalize_abcabc_input(certificate.input))) {
        ABCABCResult r(certificate);
        r.decision=Decision::Indeterminate;
        r.termination=TerminationReason::BackendFailure;
        return r;
    }
    // Do not trust the stored decision/reason; recompute theorem conditions from input.
    return result_from_eval(certificate.input,evaluate(certificate.input));
}

} // namespace pcg
