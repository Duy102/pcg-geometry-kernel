#include <pcg/core.hpp>
#include <pcg/geometry.hpp>
#include <pcg/abcabc.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
int failures = 0;
void check(bool cond, const std::string& name) {
    if (!cond) { std::cerr << "FAIL: " << name << "\n"; ++failures; }
}

template<class F> void check_domain_error(F&& f, pcg::DomainError expected, const std::string& name) {
    try { (void)f(); check(false, name + " (no error)"); }
    catch (const pcg::DomainException& e) { check(e.code() == expected, name); }
    catch (...) { check(false, name + " (wrong exception)"); }
}
}

int main() {
    using pcg::PiRational;

    // Phase 1: exact pi-multiple and domain semantics.
    check(pcg::sin_pi_sign(PiRational{1, 3}) == 1, "sin(pi/3) exact sign positive");
    check(pcg::sin_pi_sign(PiRational{4, 3}) == -1, "sin(4pi/3) exact sign negative");
    check(pcg::sin_pi_sign(PiRational{2, 1}) == 0, "sin(2pi) exact zero");
    check_domain_error([] { return pcg::Turn{PiRational{0,1}}; }, pcg::DomainError::ZeroTurnOutsideTheoremDomain, "zero turn rejected");
    check_domain_error([] { return pcg::Turn{PiRational{2,1}}; }, pcg::DomainError::FullOrOverTurnOutsideTheoremDomain, "full turn rejected");

    auto s = pcg::certified_sin_pi(PiRational{1,6});
    check(s.lower() <= 0.5 && s.upper() >= 0.5, "certified sin(pi/6) encloses 0.5");
    auto c = pcg::certified_cos_pi(PiRational{1,3});
    check(c.lower() <= 0.5 && c.upper() >= 0.5, "certified cos(pi/3) encloses 0.5");

    // Phase 2: generic shared-start relation and PCG ExtraHit are distinct layers.
    pcg::SharedStartArc a{PiRational{0,1}, pcg::Turn{PiRational{1,2}}, pcg::PositiveInterval::point(1.0)};
    pcg::SharedStartArc b{PiRational{1,2}, pcg::Turn{PiRational{1,2}}, pcg::PositiveInterval::point(1.0)};
    auto generic = pcg::shared_start_second_intersection(a,b);
    check(generic.status != pcg::GeometricStatus::BackendFailure, "generic shared-start oracle executes");

    // Phase 3: theorem-backed rational witness for the all-positive r=2 orbit.
    // Source uses x_i=tau_i/(2*pi)=(1/6,1/2,1/6,1/2,1/6,1/2),
    // therefore tau_i/pi=(1/3,1,1/3,1,1/3,1).
    pcg::ABCABCInput good{{
        pcg::Turn{PiRational{1,3}}, pcg::Turn{PiRational{1,1}}, pcg::Turn{PiRational{1,3}},
        pcg::Turn{PiRational{1,1}}, pcg::Turn{PiRational{1,3}}, pcg::Turn{PiRational{1,1}}
    }};
    auto solved = pcg::solve_abcabc(good);
    check(solved.decision == pcg::Decision::Realizable, "ABCABC equal 2pi/3 is realizable");
    check(solved.assurance == pcg::ArithmeticAssurance::CertifiedNumerical || solved.assurance == pcg::ArithmeticAssurance::Exact,
          "REALIZABLE has certified assurance");
    auto verified = pcg::verify_abcabc_certificate(solved.certificate);
    check(verified.decision == pcg::Decision::Realizable, "ABCABC certificate independently verifies");
    auto canonical = pcg::canonicalize_abcabc_input(good);
    auto digest = pcg::sha256_hex(canonical);
    check(digest.size() == 64, "canonical ABCABC input has SHA-256 digest");
    check(solved.certificate.canonical_input_digest == digest, "certificate binds canonical input digest");
    auto json1 = pcg::serialize_abcabc_certificate(solved.certificate);
    auto json2 = pcg::serialize_abcabc_certificate(solved.certificate);
    check(json1 == json2, "certificate serialization deterministic");

    auto wrong_theorem = solved.certificate;
    wrong_theorem.theorem_id = "PCG-ABCABC-THM-TAMPERED";
    auto wrong_theorem_verified = pcg::verify_abcabc_certificate(wrong_theorem);
    check(wrong_theorem_verified.decision == pcg::Decision::Indeterminate &&
          wrong_theorem_verified.termination == pcg::TerminationReason::BackendFailure,
          "certificate verifier rejects tampered theorem identity");

    auto wrong_source = solved.certificate;
    wrong_source.source_digest.value = std::string(64, '0');
    auto wrong_source_verified = pcg::verify_abcabc_certificate(wrong_source);
    check(wrong_source_verified.decision == pcg::Decision::Indeterminate &&
          wrong_source_verified.termination == pcg::TerminationReason::BackendFailure,
          "certificate verifier rejects tampered theorem source digest");

    // Symmetric all-2pi/3 input is numerically/geometrically degenerate for the
    // current certified cross-passage path and must not be overclaimed.
    pcg::ABCABCInput symmetric{{
        pcg::Turn{PiRational{2,3}}, pcg::Turn{PiRational{2,3}}, pcg::Turn{PiRational{2,3}},
        pcg::Turn{PiRational{2,3}}, pcg::Turn{PiRational{2,3}}, pcg::Turn{PiRational{2,3}}
    }};
    auto sym = pcg::solve_abcabc(symmetric);
    check(sym.decision == pcg::Decision::Indeterminate, "degenerate symmetric ABCABC is not overclaimed");

    // Exact phase obstruction.
    pcg::ABCABCInput bad_phase = good;
    bad_phase.turns[0] = pcg::Turn{PiRational{1,2}};
    auto rejected = pcg::solve_abcabc(bad_phase);
    check(rejected.decision == pcg::Decision::NotRealizable, "ABCABC phase failure rejected");
    check(rejected.assurance == pcg::ArithmeticAssurance::Exact, "phase obstruction is exact");
    check(pcg::verify_abcabc_certificate(rejected.certificate).decision == pcg::Decision::NotRealizable,
          "phase obstruction certificate verifies");

    if (failures) {
        std::cerr << failures << " test(s) failed\n";
        return 1;
    }
    std::cout << "All PCG tests passed\n";
    return 0;
}
