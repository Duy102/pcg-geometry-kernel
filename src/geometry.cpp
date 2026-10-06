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


namespace {

CertifiedTruth closed_nonnegative(const I& x) {
    if (x.lower() >= 0.0) return CertifiedTruth::True;
    if (x.upper() < 0.0) return CertifiedTruth::False;
    return CertifiedTruth::Indeterminate;
}

CertifiedTruth truth_and(CertifiedTruth a, CertifiedTruth b) {
    if (a==CertifiedTruth::False || b==CertifiedTruth::False)
        return CertifiedTruth::False;
    if (a==CertifiedTruth::True && b==CertifiedTruth::True)
        return CertifiedTruth::True;
    return CertifiedTruth::Indeterminate;
}

CertifiedTruth truth_or(CertifiedTruth a, CertifiedTruth b) {
    if (a==CertifiedTruth::True || b==CertifiedTruth::True)
        return CertifiedTruth::True;
    if (a==CertifiedTruth::False && b==CertifiedTruth::False)
        return CertifiedTruth::False;
    return CertifiedTruth::Indeterminate;
}

I cross_i(const I& ax,const I& ay,const I& bx,const I& by) {
    return ax*by-ay*bx;
}

struct FiniteCenterRadius {
    I cx;
    I cy;
    I radius;
};

FiniteCenterRadius finite_center_radius(const CertifiedFiniteArc& arc) {
    const I sin_half=certified_sin_pi(arc.turn.pi/2).interval();
    const I rho=arc.chord.scalar().interval()/(I(2.0)*sin_half);
    const I s=certified_sin_pi(arc.tangent_phase_pi).interval();
    const I c=certified_cos_pi(arc.tangent_phase_pi).interval();
    const I cx=arc.source.x.interval()-rho*s;
    const I cy=arc.source.y.interval()+rho*c;
    return {cx,cy,abs_i(rho)};
}

CertifiedTruth finite_arc_membership(const CertifiedFiniteArc& arc,
                                     const FiniteCenterRadius& support,
                                     const I& qx,const I& qy) {
    const I rmx=arc.source.x.interval()-support.cx;
    const I rmy=arc.source.y.interval()-support.cy;
    const I rpx=arc.target.x.interval()-support.cx;
    const I rpy=arc.target.y.interval()-support.cy;
    const I xx=qx-support.cx;
    const I xy=qy-support.cy;
    const double orient=arc.turn.pi.value > Rational(BigInt{0}) ? 1.0 : -1.0;

    const I A=I(orient)*cross_i(rmx,rmy,xx,xy);
    const I B=I(orient)*cross_i(xx,xy,rpx,rpy);
    const auto a=closed_nonnegative(A);
    const auto b=closed_nonnegative(B);

    Rational mag=arc.turn.pi.value;
    if (mag<Rational(BigInt{0})) mag=-mag;
    if (mag<Rational(BigInt{1})) return truth_and(a,b);
    if (mag==Rational(BigInt{1})) return a;
    return truth_or(a,b);
}

CertifiedTruth common_membership(const CertifiedFiniteArc& a,
                                 const FiniteCenterRadius& ca,
                                 const CertifiedFiniteArc& b,
                                 const FiniteCenterRadius& cb,
                                 const I& qx,const I& qy) {
    return truth_and(finite_arc_membership(a,ca,qx,qy),
                     finite_arc_membership(b,cb,qx,qy));
}

} // namespace

CertifiedTruth pcg_hit(const CertifiedFiniteArc& a, const CertifiedFiniteArc& b) {
    try {
        const auto ca=finite_center_radius(a);
        const auto cb=finite_center_radius(b);

        const I dx=cb.cx-ca.cx;
        const I dy=cb.cy-ca.cy;
        const I d2=dx*dx+dy*dy;

        if (d2.upper()<=0.0) {
            const I dr=abs_i(ca.radius-cb.radius);
            if (dr.lower()>0.0) return CertifiedTruth::False;
            return CertifiedTruth::Indeterminate;
        }
        if (d2.lower()<=0.0) return CertifiedTruth::Indeterminate;

        const I d=boost::numeric::sqrt(d2);
        const I sum=ca.radius+cb.radius;
        const I diff=abs_i(ca.radius-cb.radius);

        if (d.lower()>sum.upper()) return CertifiedTruth::False;
        if (d.upper()<diff.lower()) return CertifiedTruth::False;

        // Phase 6A deliberately requires a certified two-intersection support
        // configuration. Tangency boundaries are left Indeterminate.
        if (!(d.upper()<sum.lower() && d.lower()>diff.upper()))
            return CertifiedTruth::Indeterminate;

        const I along=(ca.radius*ca.radius-cb.radius*cb.radius+d2)/(I(2.0)*d);
        const I h2=ca.radius*ca.radius-along*along;
        if (h2.lower()<=0.0) return CertifiedTruth::Indeterminate;
        const I h=boost::numeric::sqrt(h2);

        const I ux=dx/d;
        const I uy=dy/d;
        const I bx=ca.cx+along*ux;
        const I by=ca.cy+along*uy;
        const I px=-uy;
        const I py=ux;

        const auto q1=common_membership(
            a,ca,b,cb,bx+h*px,by+h*py);
        if (q1==CertifiedTruth::True) return CertifiedTruth::True;

        const auto q2=common_membership(
            a,ca,b,cb,bx-h*px,by-h*py);
        if (q2==CertifiedTruth::True) return CertifiedTruth::True;

        if (q1==CertifiedTruth::False && q2==CertifiedTruth::False)
            return CertifiedTruth::False;
        return CertifiedTruth::Indeterminate;
    } catch (...) {
        return CertifiedTruth::Indeterminate;
    }
}

} // namespace pcg
