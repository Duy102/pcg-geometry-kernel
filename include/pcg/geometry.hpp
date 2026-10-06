#pragma once
#include <pcg/core.hpp>
#include <optional>

namespace pcg {

struct SupportingCircle {
    CertifiedScalar cx;
    CertifiedScalar cy;
    CertifiedScalar radius;
};

struct CircularArc {
    Point2 source;
    Point2 target;
    Turn turn;
    SupportingCircle support;
};

enum class CertifiedTruth { True, False, Indeterminate };
enum class GeometricStatus {
    RegularSecondIntersection,
    TangentAtStart,
    CoincidentSupport,
    Indeterminate,
    BackendFailure
};

enum class SupportRelation {
    Distinct,
    Coincident,
    Indeterminate,
    BackendFailure
};

struct SharedStartArc {
    PiRational tangent_phase_pi;
    Turn turn;
    PositiveInterval chord;
};

struct SharedStartIntersection {
    GeometricStatus status{GeometricStatus::Indeterminate};
    CertifiedScalar qx;
    CertifiedScalar qy;
};

SupportRelation shared_start_support_relation(const SharedStartArc& a, const SharedStartArc& b);
SharedStartIntersection shared_start_second_intersection(const SharedStartArc& a, const SharedStartArc& b);
CertifiedTruth shared_start_arc_contains_second(const SharedStartArc& arc,
                                                PiRational normalized_tangent_phase,
                                                const CertifiedScalar& qx,
                                                const CertifiedScalar& qy);
CertifiedTruth pcg_extra_hit(const SharedStartArc& a, const SharedStartArc& b);

} // namespace pcg
