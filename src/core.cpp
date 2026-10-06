#include <pcg/core.hpp>
#include <boost/numeric/interval/transc.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <limits>
#include <openssl/sha.h>
#include <iomanip>
#include <sstream>
#include <utility>

namespace pcg {
namespace {
using I = detail::Interval;

I pi_interval() {
    // Adjacent binary64 values rigorously bracketing mathematical pi.
    return I(0x1.921fb54442d17p+1, 0x1.921fb54442d19p+1);
}

I rational_interval(const Rational& q) {
    // Convert an arbitrary-precision exact rational through a high-precision
    // decimal approximation, then widen by one binary64 ulp on each side.
    // The exact rational remains the theorem-facing representation; this
    // conversion is used only to seed the certified interval substrate.
    using Dec100 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>>;
    const Dec100 approx = Dec100(q.numerator()) / Dec100(q.denominator());
    const double d = approx.convert_to<double>();
    return I(std::nextafter(d, -std::numeric_limits<double>::infinity()),
             std::nextafter(d,  std::numeric_limits<double>::infinity()));
}

BigInt floor_rational(const Rational& r) {
    const BigInt& n = r.numerator();
    const BigInt& d = r.denominator();
    BigInt q = n / d;
    const BigInt rem = n % d;
    if (rem != 0 && n < 0) --q;
    return q;
}
}

PiRational::PiRational(std::int64_t n, std::int64_t d)
    : PiRational(BigInt{n}, BigInt{d}) {}

PiRational::PiRational(BigInt n, BigInt d)
    : value(std::move(n), std::move(d)) {}

PiRational::PiRational(Rational r)
    : value(std::move(r)) {}

PiRational operator+(PiRational a, PiRational b) { return PiRational{a.value+b.value}; }
PiRational operator-(PiRational a, PiRational b) { return PiRational{a.value-b.value}; }
PiRational operator-(PiRational a) { return PiRational{-a.value}; }
PiRational operator*(PiRational a, Rational b) { return PiRational{a.value*b}; }
PiRational operator/(PiRational a, std::int64_t d) { return PiRational{a.value/Rational(BigInt{d})}; }
bool operator==(PiRational a, PiRational b) { return a.value==b.value; }

int sin_pi_sign(PiRational x) {
    // Reduce q modulo 2 into [0,2), exactly.
    Rational two(BigInt{2});
    const auto k = floor_rational(x.value / two);
    Rational r = x.value - Rational(BigInt{2}*k);
    if (r == Rational(BigInt{0}) || r == Rational(BigInt{1})) return 0;
    return r < Rational(BigInt{1}) ? 1 : -1;
}

bool is_even_integer(PiRational x) {
    return x.value.denominator() == 1 && (x.value.numerator() % 2 == 0);
}

CertifiedScalar::CertifiedScalar() : value_(0.0) {}
CertifiedScalar::CertifiedScalar(double lo, double hi) : value_(lo,hi) {}
CertifiedScalar::CertifiedScalar(detail::Interval interval) : value_(interval) {}
double CertifiedScalar::lower() const { return value_.lower(); }
double CertifiedScalar::upper() const { return value_.upper(); }
const detail::Interval& CertifiedScalar::interval() const { return value_; }

PositiveInterval PositiveInterval::point(double x) {
    if (!(x > 0.0) || !std::isfinite(x)) throw std::invalid_argument("positive finite value required");
    return PositiveInterval(CertifiedScalar{x,x});
}
PositiveInterval::PositiveInterval(CertifiedScalar v) : value_(std::move(v)) {
    if (!(value_.lower() > 0.0)) throw std::invalid_argument("interval must be strictly positive");
}
const CertifiedScalar& PositiveInterval::scalar() const { return value_; }

CertifiedScalar certified_radians(PiRational x) {
    return CertifiedScalar(rational_interval(x.value) * pi_interval());
}
CertifiedScalar certified_sin_pi(PiRational x) {
    return CertifiedScalar(boost::numeric::sin(certified_radians(x).interval()));
}
CertifiedScalar certified_cos_pi(PiRational x) {
    return CertifiedScalar(boost::numeric::cos(certified_radians(x).interval()));
}
CertifiedScalar certified_sinc(PiRational x) {
    if (x.value == Rational(0)) return CertifiedScalar(1.0,1.0);
    auto xr = certified_radians(x).interval();
    return CertifiedScalar(boost::numeric::sin(xr) / xr);
}
CertifiedScalar certified_abs(const CertifiedScalar& x) {
    const auto lo=x.lower(), hi=x.upper();
    if (lo >= 0.0) return x;
    if (hi <= 0.0) return CertifiedScalar(-hi,-lo);
    return CertifiedScalar(0.0,std::max(-lo,hi));
}

std::string rational_string(const PiRational& q) {
    std::ostringstream os;
    os << q.value.numerator() << "/" << q.value.denominator();
    return os.str();
}

std::string sha256_hex(const std::string& input) {
    unsigned char digest[SHA256_DIGEST_LENGTH];
    SHA256(reinterpret_cast<const unsigned char*>(input.data()), input.size(), digest);
    std::ostringstream os;
    os << std::hex << std::setfill('0');
    for (unsigned char b : digest) os << std::setw(2) << static_cast<unsigned int>(b);
    return os.str();
}

DomainException::DomainException(DomainError code, const std::string& what)
    : std::runtime_error(what), code_(code) {}
DomainError DomainException::code() const noexcept { return code_; }

Turn::Turn(PiRational q) : pi(q) {
    if (q.value == Rational(0))
        throw DomainException(DomainError::ZeroTurnOutsideTheoremDomain, "turn must be nonzero");
    Rational a = q.value < Rational(0) ? -q.value : q.value;
    if (a >= Rational(2))
        throw DomainException(DomainError::FullOrOverTurnOutsideTheoremDomain, "turn must satisfy |tau| < 2pi");
}

} // namespace pcg
