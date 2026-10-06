#include <pcg/core.hpp>
#include <boost/numeric/interval/transc.hpp>
#include <boost/numeric/interval/utility.hpp>
#include <cmath>
#include <limits>
#include <openssl/sha.h>
#include <iomanip>
#include <sstream>

namespace pcg {
namespace {
using I = detail::Interval;

I pi_interval() {
    // Adjacent binary64 values rigorously bracketing mathematical pi.
    return I(0x1.921fb54442d17p+1, 0x1.921fb54442d19p+1);
}

I rational_interval(const Rational& q) {
    const long double exactish = static_cast<long double>(q.numerator()) / static_cast<long double>(q.denominator());
    const double d = static_cast<double>(exactish);
    return I(std::nextafter(d, -std::numeric_limits<double>::infinity()),
             std::nextafter(d,  std::numeric_limits<double>::infinity()));
}

std::int64_t floor_rational(const Rational& r) {
    auto n = r.numerator();
    auto d = r.denominator();
    auto q = n / d;
    auto rem = n % d;
    if (rem != 0 && n < 0) --q;
    return q;
}
}

PiRational::PiRational(std::int64_t n, std::int64_t d) : value(n,d) {
    if (d == 0) throw std::invalid_argument("zero denominator");
}
PiRational operator+(PiRational a, PiRational b) { return PiRational{(a.value+b.value).numerator(), (a.value+b.value).denominator()}; }
PiRational operator-(PiRational a, PiRational b) { auto r=a.value-b.value; return PiRational{r.numerator(),r.denominator()}; }
PiRational operator-(PiRational a) { return PiRational{-a.value.numerator(),a.value.denominator()}; }
PiRational operator*(PiRational a, Rational b) { auto r=a.value*b; return PiRational{r.numerator(),r.denominator()}; }
PiRational operator/(PiRational a, std::int64_t d) { auto r=a.value/Rational(d); return PiRational{r.numerator(),r.denominator()}; }
bool operator==(PiRational a, PiRational b) { return a.value==b.value; }

int sin_pi_sign(PiRational x) {
    // Reduce q modulo 2 into [0,2), exactly.
    Rational two(2);
    const auto k = floor_rational(x.value / two);
    Rational r = x.value - Rational(2*k);
    if (r == Rational(0) || r == Rational(1)) return 0;
    return r < Rational(1) ? 1 : -1;
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
    return std::to_string(q.value.numerator()) + "/" + std::to_string(q.value.denominator());
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
