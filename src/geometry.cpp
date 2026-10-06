#include <pcg/geometry.hpp>
#include <boost/numeric/interval/transc.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <algorithm>
#include <cmath>
#include <limits>

namespace pcg {
namespace {
using I = detail::Interval;

CertifiedTruth strict_positive(const I& x) {
    if (x.lower() > 0.0) return CertifiedTruth::True;
    if (x.upper() <= 0.0) return CertifiedTruth::False;
    return CertifiedTruth::Indeterminate;
}

CertifiedTruth strict_less(const I& a, const I& b) {
    if (a.upper() < b.lower()) return CertifiedTruth::True;
    if (a.lower() >= b.upper()) return CertifiedTruth::False;
    return CertifiedTruth::Indeterminate;
}

I abs_i(const I& x) {
    if (x.lower() >= 0.0) return x;
    if (x.upper() <= 0.0) return I(-x.upper(), -x.lower());
    return I(0.0, std::max(-x.lower(), x.upper()));
}

std::optional<I> atan2_i(const I& y, const I& x) {
    const I pi = certified_radians(PiRational{1,1}).interval();
    if (x.lower() > 0.0) {
        return boost::numeric::atan(y / x);
    }
    if (x.upper() < 0.0) {
        // The negative x-axis is the principal-argument branch cut.
        // If y contains zero, including the exact antipodal case, do not
        // silently choose +pi or -pi: theorem-facing callers must return
        // INDETERMINATE unless a separate domain certificate resolves it.
        if (y.lower() <= 0.0 && y.upper() >= 0.0) return std::nullopt;
        if (y.lower() > 0.0) return boost::numeric::atan(y / x) + pi;
        if (y.upper() < 0.0) return boost::numeric::atan(y / x) - pi;
        return std::nullopt;
    }
    return std::nullopt;
}

struct CenterAndRadius {
    I cx; I cy; I signed_radius;
};

CenterAndRadius center_for_normalized_arc(const SharedStartArc& arc, PiRational relative_phase) {
    const I sin_half = certified_sin_pi(arc.turn.pi / 2).interval();
    const I rho = arc.chord.scalar().interval() / (I(2.0) * sin_half);
    const I s = certified_sin_pi(relative_phase).interval();
    const I c = certified_cos_pi(relative_phase).interval();
    return { -rho*s, rho*c, rho };
}

CertifiedTruth membership_from_q(const SharedStartArc& arc,
                                 PiRational relative_phase,
                                 const I& qx,
                                 const I& qy) {
    // Rotate Q into the arc's initial-tangent frame.
    const I s = certified_sin_pi(relative_phase).interval();
    const I c = certified_cos_pi(relative_phase).interval();
    const I xr = qx*c + qy*s;
    const I yr = -qx*s + qy*c;
    const auto arg = atan2_i(yr,xr);
    if (!arg) return CertifiedTruth::Indeterminate;
    const I delta = I(2.0) * *arg;
    const I tau = certified_radians(arc.turn.pi).interval();

    const auto same_sign = strict_positive(delta*tau);
    const auto within = strict_less(abs_i(delta), abs_i(tau));
    if (same_sign == CertifiedTruth::False || within == CertifiedTruth::False) return CertifiedTruth::False;
    if (same_sign == CertifiedTruth::True && within == CertifiedTruth::True) return CertifiedTruth::True;
    return CertifiedTruth::Indeterminate;
}
}

SharedStartIntersection shared_start_second_intersection(const SharedStartArc& a, const SharedStartArc& b) {
    try {
        const PiRational delta = b.tangent_phase_pi - a.tangent_phase_pi;
        const auto ca = center_for_normalized_arc(a, PiRational{0,1});
        const auto cb = center_for_normalized_arc(b, delta);
        const I wx = cb.cx - ca.cx;
        const I wy = cb.cy - ca.cy;
        const I norm2 = wx*wx + wy*wy;
        if (norm2.upper() <= 0.0) return {GeometricStatus::CoincidentSupport, {}, {}};
        if (norm2.lower() <= 0.0) return {GeometricStatus::Indeterminate, {}, {}};
        const I norm = boost::numeric::sqrt(norm2);
        const I ux = -wy / norm;
        const I uy =  wx / norm;
        const I projection = ca.cx*ux + ca.cy*uy;
        const I scale = I(2.0) * projection;
        const I qx = scale*ux;
        const I qy = scale*uy;
        const I qnorm2 = qx*qx + qy*qy;
        if (qnorm2.upper() <= 0.0) return {GeometricStatus::TangentAtStart, CertifiedScalar(qx), CertifiedScalar(qy)};
        if (qnorm2.lower() <= 0.0) return {GeometricStatus::Indeterminate, CertifiedScalar(qx), CertifiedScalar(qy)};
        return {GeometricStatus::RegularSecondIntersection, CertifiedScalar(qx), CertifiedScalar(qy)};
    } catch (...) {
        return {GeometricStatus::BackendFailure, {}, {}};
    }
}

CertifiedTruth shared_start_arc_contains_second(const SharedStartArc& arc,
                                                PiRational normalized_tangent_phase,
                                                const CertifiedScalar& qx,
                                                const CertifiedScalar& qy) {
    try {
        return membership_from_q(arc, normalized_tangent_phase, qx.interval(), qy.interval());
    } catch (...) {
        return CertifiedTruth::Indeterminate;
    }
}

CertifiedTruth pcg_extra_hit(const SharedStartArc& a, const SharedStartArc& b) {
    const auto inter = shared_start_second_intersection(a,b);
    if (inter.status == GeometricStatus::TangentAtStart) return CertifiedTruth::False;
    if (inter.status != GeometricStatus::RegularSecondIntersection) return CertifiedTruth::Indeterminate;
    const PiRational delta = b.tangent_phase_pi - a.tangent_phase_pi;
    const auto ina = membership_from_q(a, PiRational{0,1}, inter.qx.interval(), inter.qy.interval());
    const auto inb = membership_from_q(b, delta, inter.qx.interval(), inter.qy.interval());
    if (ina == CertifiedTruth::False || inb == CertifiedTruth::False) return CertifiedTruth::False;
    if (ina == CertifiedTruth::True && inb == CertifiedTruth::True) return CertifiedTruth::True;
    return CertifiedTruth::Indeterminate;
}

SupportRelation shared_start_support_relation(const SharedStartArc& a, const SharedStartArc& b) {
    try {
        const PiRational delta = b.tangent_phase_pi - a.tangent_phase_pi;

        // Exact fast path for the common same-parameter coincidence that
        // interval dependency would otherwise widen around zero.
        const auto& ac = a.chord.scalar();
        const auto& bc = b.chord.scalar();
        if (is_even_integer(delta) &&
            a.turn.pi == b.turn.pi &&
            ac.lower() == bc.lower() &&
            ac.upper() == bc.upper()) {
            return SupportRelation::Coincident;
        }

        const auto ca = center_for_normalized_arc(a, PiRational{0,1});
        const auto cb = center_for_normalized_arc(b, delta);
        const I wx = cb.cx - ca.cx;
        const I wy = cb.cy - ca.cy;
        const I norm2 = wx*wx + wy*wy;
        if (norm2.upper() <= 0.0) return SupportRelation::Coincident;
        if (norm2.lower() > 0.0) return SupportRelation::Distinct;
        return SupportRelation::Indeterminate;
    } catch (...) {
        return SupportRelation::BackendFailure;
    }
}

} // namespace pcg
