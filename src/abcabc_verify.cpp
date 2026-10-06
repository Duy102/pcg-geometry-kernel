#include <pcg/abcabc.hpp>
#include <array>
#include <utility>

namespace pcg {
namespace {

// This translation unit intentionally duplicates the top-level theorem
// recomputation instead of calling solve_abcabc() or its private evaluate()
// path. It still shares the documented geometry/numerical TCB.

PiRational vhalf(PiRational q) { return q / 2; }

std::array<PiRational,6> vqturns(const ABCABCInput& in) {
    std::array<PiRational,6> q{};
    for (std::size_t i=0;i<6;++i) q[i]=in.turns[i].pi;
    return q;
}

bool vphase_ok(const std::array<PiRational,6>& q) {
    const PiRational p1 = vhalf(q[0]) + q[1] + q[2] + vhalf(q[3]);
    const PiRational p2 = vhalf(q[1]) + q[2] + q[3] + vhalf(q[4]);
    const PiRational p3 = vhalf(q[2]) + q[3] + q[4] + vhalf(q[5]);
    return is_even_integer(p1) && is_even_integer(p2) && is_even_integer(p3);
}

std::array<PiRational,6> vphases(const std::array<PiRational,6>& q) {
    std::array<PiRational,6> p{};
    p[0]=PiRational{0,1};
    for (std::size_t i=1;i<6;++i) p[i]=p[i-1]+q[i-1];
    return p;
}

std::array<PiRational,3> valphas(const std::array<PiRational,6>& q,
                                 const std::array<PiRational,6>& p) {
    return {p[0]+vhalf(q[0]), p[1]+vhalf(q[1]), p[2]+vhalf(q[2])};
}

bool vpositive_closure_ok(const std::array<PiRational,3>& a) {
    const int s1=sin_pi_sign(a[2]-a[1]);
    const int s2=sin_pi_sign(a[0]-a[2]);
    const int s3=sin_pi_sign(a[1]-a[0]);
    return s1!=0 && s1==s2 && s2==s3;
}

PositiveInterval vchord(PiRational q) {
    return PositiveInterval(certified_abs(certified_sin_pi(q)));
}

SharedStartArc vforward(std::size_t i,
                        const std::array<PiRational,6>& p,
                        const ABCABCInput& in,
                        const std::array<PositiveInterval,6>& c) {
    return SharedStartArc{p[i],in.turns[i],c[i]};
}

SharedStartArc vreverse(std::size_t i,
                        const std::array<PiRational,6>& p,
                        const ABCABCInput& in,
                        const std::array<PositiveInterval,6>& c) {
    return SharedStartArc{p[i]+in.turns[i].pi+PiRational{1,1},
                          Turn{-in.turns[i].pi},c[i]};
}

using VPair = std::pair<SharedStartArc,SharedStartArc>;

std::array<VPair,6> vcross_pairs(const std::array<PiRational,6>& p,
                                 const ABCABCInput& in,
                                 const std::array<PositiveInterval,6>& c) {
    return {{
        {vforward(0,p,in,c), vreverse(2,p,in,c)},
        {vforward(3,p,in,c), vreverse(5,p,in,c)},
        {vforward(1,p,in,c), vreverse(3,p,in,c)},
        {vforward(4,p,in,c), vreverse(0,p,in,c)},
        {vforward(2,p,in,c), vreverse(4,p,in,c)},
        {vforward(5,p,in,c), vreverse(1,p,in,c)}
    }};
}

struct VEval {
    Decision decision;
    ProofKind kind;
    ArithmeticAssurance assurance;
    ABCABCProofReason reason;
};

VEval verify_evaluate(const ABCABCInput& in) {
    const auto q=vqturns(in);
    if (!vphase_ok(q))
        return {Decision::NotRealizable,ProofKind::TheoremInstantiation,
                ArithmeticAssurance::Exact,ABCABCProofReason::PhaseCongruenceObstruction};

    const auto p=vphases(q);
    const auto a=valphas(q,p);
    if (!vpositive_closure_ok(a))
        return {Decision::NotRealizable,ProofKind::TheoremInstantiation,
                ArithmeticAssurance::Exact,ABCABCProofReason::PositiveClosureObstruction};

    const auto genericity=certify_abcabc_genericity(in);
    if (genericity==ABCABCGenericityStatus::Violated)
        return {Decision::Unsupported,ProofKind::FinitePredicateCertificate,
                ArithmeticAssurance::CertifiedNumerical,ABCABCProofReason::TheoremDomainViolation};
    if (genericity!=ABCABCGenericityStatus::Satisfied)
        return {Decision::Indeterminate,ProofKind::None,
                ArithmeticAssurance::None,ABCABCProofReason::NumericalIndeterminacy};

    std::array<PositiveInterval,6> c{
        vchord(a[2]-a[1]),
        vchord(a[0]-a[2]),
        vchord(a[1]-a[0]),
        vchord(a[2]-a[1]),
        vchord(a[0]-a[2]),
        vchord(a[1]-a[0])
    };

    bool indeterminate=false;
    for (const auto& [x,y] : vcross_pairs(p,in,c)) {
        const auto inter=shared_start_second_intersection(x,y);
        if (inter.status!=GeometricStatus::RegularSecondIntersection) {
            indeterminate=true;
            continue;
        }

        const PiRational delta=y.tangent_phase_pi-x.tangent_phase_pi;
        const auto in_x=shared_start_arc_contains_second(
            x,PiRational{0,1},inter.qx,inter.qy);
        const auto in_y=shared_start_arc_contains_second(
            y,delta,inter.qx,inter.qy);

        if (in_x==CertifiedTruth::False || in_y==CertifiedTruth::False)
            continue;

        if (in_x==CertifiedTruth::True && in_y==CertifiedTruth::True)
            return {Decision::NotRealizable,ProofKind::FinitePredicateCertificate,
                    ArithmeticAssurance::CertifiedNumerical,ABCABCProofReason::ExtraHitObstruction};

        indeterminate=true;
    }

    if (indeterminate)
        return {Decision::Indeterminate,ProofKind::None,
                ArithmeticAssurance::None,ABCABCProofReason::NumericalIndeterminacy};

    return {Decision::Realizable,ProofKind::CompositeCertificate,
            ArithmeticAssurance::CertifiedNumerical,ABCABCProofReason::AllConditionsSatisfied};
}

ABCABCResult invalid_certificate(const ABCABCCertificate& c) {
    ABCABCResult r(c);
    r.decision=Decision::Indeterminate;
    r.proof_kind=ProofKind::None;
    r.assurance=ArithmeticAssurance::None;
    r.termination=TerminationReason::BackendFailure;
    return r;
}

ABCABCResult result_from_reference(const ABCABCCertificate& c, const VEval& e) {
    return ABCABCResult{
        e.decision,e.kind,e.assurance,
        e.decision==Decision::Indeterminate
            ? TerminationReason::PrecisionLimit
            : TerminationReason::Completed,
        c
    };
}

} // namespace

ABCABCResult verify_abcabc_certificate(const ABCABCCertificate& certificate) {
    if (certificate.schema!="pcg-abcabc-certificate" ||
        certificate.schema_version!="1.0" ||
        certificate.theorem_id!="PCG-ABCABC-THM-001" ||
        certificate.source_digest.algorithm!="SHA-256" ||
        certificate.source_digest.value!="2aa832a00031035d964c0d9362e4409861e19b97a9d764f3123a42edef400beb" ||
        certificate.canonical_input_digest!=sha256_hex(canonicalize_abcabc_input(certificate.input))) {
        return invalid_certificate(certificate);
    }

    VEval e;
    try {
        e=verify_evaluate(certificate.input);
    } catch (...) {
        return invalid_certificate(certificate);
    }

    // The certificate is evidence, not merely an input container. Tampering
    // with its claimed reason/proof/assurance must be detectable even when the
    // bound input itself is unchanged.
    if (certificate.reason!=e.reason ||
        certificate.proof_kind!=e.kind ||
        certificate.assurance!=e.assurance) {
        return invalid_certificate(certificate);
    }

    return result_from_reference(certificate,e);
}

} // namespace pcg
