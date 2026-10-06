#pragma once

#include <boost/numeric/interval.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>
#include <array>
#include <cstdint>
#include <stdexcept>
#include <string>

namespace pcg {

using BigInt = boost::multiprecision::cpp_int;
using Rational = boost::rational<BigInt>;

struct PiRational {
    Rational value;
    PiRational(std::int64_t n = 0, std::int64_t d = 1);
    PiRational(BigInt n, BigInt d);
    explicit PiRational(Rational r);
};

PiRational operator+(PiRational a, PiRational b);
PiRational operator-(PiRational a, PiRational b);
PiRational operator-(PiRational a);
PiRational operator*(PiRational a, Rational b);
PiRational operator/(PiRational a, std::int64_t d);
bool operator==(PiRational a, PiRational b);

int sin_pi_sign(PiRational x);
bool is_even_integer(PiRational x);

namespace detail {
using Rounding = boost::numeric::interval_lib::save_state<
    boost::numeric::interval_lib::rounded_transc_std<double,
        boost::numeric::interval_lib::rounded_arith_std<double>>>;
using Policies = boost::numeric::interval_lib::policies<
    Rounding,
    boost::numeric::interval_lib::checking_strict<double>>;
using Interval = boost::numeric::interval<double, Policies>;
}

class CertifiedScalar {
public:
    CertifiedScalar();
    CertifiedScalar(double lo, double hi);
    explicit CertifiedScalar(detail::Interval interval);
    double lower() const;
    double upper() const;
    const detail::Interval& interval() const;
private:
    detail::Interval value_;
};

class PositiveInterval {
public:
    static PositiveInterval point(double x);
    explicit PositiveInterval(CertifiedScalar v);
    const CertifiedScalar& scalar() const;
private:
    CertifiedScalar value_;
};

CertifiedScalar certified_radians(PiRational x);
CertifiedScalar certified_sin_pi(PiRational x);
CertifiedScalar certified_cos_pi(PiRational x);
CertifiedScalar certified_sinc(PiRational radians_over_pi);
CertifiedScalar certified_abs(const CertifiedScalar& x);

struct Angle { PiRational pi; };

struct Point2 { double x{}; double y{}; };
struct Vector2 { double x{}; double y{}; };

struct CanonicalFrame {
    Point2 origin{0.0, 0.0};
    PiRational rotation{0,1};
    double scale{1.0};
};

enum class InputError {
    InvalidSchema,
    NonFiniteInput,
    InvalidIndex,
    UnsupportedSchemaVersion
};

enum class DomainError {
    ZeroTurnOutsideTheoremDomain,
    FullOrOverTurnOutsideTheoremDomain,
    WrongTraceFamily,
    RequiredNondegeneracyViolated
};

class DomainException : public std::runtime_error {
public:
    DomainException(DomainError code, const std::string& what);
    DomainError code() const noexcept;
private:
    DomainError code_;
};

struct Turn {
    PiRational pi;
    explicit Turn(PiRational q);
};

enum class Decision { Realizable, NotRealizable, Indeterminate, Unsupported };
enum class ProofKind {
    TheoremInstantiation,
    ConstructiveWitness,
    FinitePredicateCertificate,
    CompositeCertificate,
    None
};
enum class ArithmeticAssurance { Exact, CertifiedNumerical, None };
enum class TerminationReason { Completed, PrecisionLimit, BackendFailure };

struct TheoremId { std::string value; };
struct SourceDigest { std::string algorithm; std::string value; };

std::string sha256_hex(const std::string& input);
std::string rational_string(const PiRational& q);

} // namespace pcg
