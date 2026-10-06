#include <pcg/qe_backend.hpp>

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

int failures=0;

void check(bool cond,const std::string& name) {
    if (!cond) {
        std::cerr << "FAIL: " << name << "\n";
        ++failures;
    }
}

std::size_t count_substring(
    const std::string& s,
    const std::string& needle) {
    std::size_t count=0;
    std::size_t pos=0;
    while ((pos=s.find(needle,pos))!=std::string::npos) {
        ++count;
        pos+=needle.size();
    }
    return count;
}

bool interval_contains(
    const pcg::TrigIsolationDefinition& d,
    const pcg::Rational& x) {
    return d.cosine_lower<x && x<d.cosine_upper;
}

}

int main() {
    using namespace pcg;

    NetworkClosureInput triangle{
        3,{0,1,2},
        {
            Turn{PiRational{1,2}},
            Turn{PiRational{1,1}},
            Turn{PiRational{3,2}}
        }
    };

    const auto source=compile_exact_semialgebraic_program(triangle);
    auto p=lower_to_qepcad_rcf_program(source);

    check(
        p.source_program_digest==
        exact_semialgebraic_program_digest(source),
        "Phase 6B2C binds exactly to the Phase 6B2B polynomial AST");

    check(
        verify_qepcad_rcf_program(p).status==
        QepcadRcfValidationStatus::Valid,
        "QEPCAD rational RCF program verifies by deterministic recompilation");

    check(
        !p.coefficient_variables.empty() &&
        p.variables.size()>source.variables.size(),
        "algebraic coefficients are lifted to quantified real variables");

    bool saw_half=false;
    bool saw_zero=false;
    for (const auto& d:p.trig_definitions) {
        if (d.phase_pi==PiRational{1,3}) {
            saw_half=true;
            check(
                interval_contains(
                    d,Rational(BigInt{1},BigInt{2})),
                "cos(pi/3)=1/2 lies in its exact Sturm isolation interval");
            check(
                d.sine_sign==1,
                "sin(pi/3) exact sign is positive");
        }
        if (d.phase_pi==PiRational{1,2}) {
            saw_zero=true;
            check(
                interval_contains(d,Rational(BigInt{0})),
                "cos(pi/2)=0 lies in its exact Sturm isolation interval");
            check(
                d.sine_sign==1,
                "sin(pi/2) exact sign is positive");
        }
    }
    check(saw_half,"adapter emitted a trig definition for phase pi/3");
    check(saw_zero,"adapter emitted a trig definition for phase pi/2");

    check(
        p.maximum_total_degree>=2 &&
        p.rational_atom_count>source.polynomial_atom_count,
        "rational RCF formula records algebraic definitions in addition to geometry atoms");

    const auto serialized=serialize_qepcad_rcf_program(p);
    check(
        serialized.find("cosine_polynomial")!=std::string::npos &&
        serialized.find("source_program_digest")!=std::string::npos,
        "canonical backend serialization includes isolation and source binding");

    const auto d1=qepcad_rcf_program_digest(p);
    const auto d2=qepcad_rcf_program_digest(
        compile_qepcad_rcf_program(triangle));
    check(
        d1==d2 && d1.size()==64,
        "QEPCAD backend package SHA-256 is deterministic");

    const auto qepcad=export_qepcad_input(p);
    check(
        qepcad.find("\n0\n")!=std::string::npos,
        "QEPCAD decision input declares zero free variables");
    check(
        qepcad.find("(E v")!=std::string::npos,
        "QEPCAD export contains existential coefficient and geometry variables");
    check(
        qepcad.find("(A v")!=std::string::npos,
        "not-exists bad-pair incidences prenex to universal QEPCAD variables");
    check(
        qepcad.find(" /\\ ")!=std::string::npos &&
        qepcad.find("finish\n")!=std::string::npos,
        "QEPCAD export uses documented Boolean syntax and finish command");
    check(
        qepcad.find("sin_pi")==std::string::npos &&
        qepcad.find("cos_pi")==std::string::npos,
        "QEPCAD payload contains only rational polynomial constraints");
    check(
        count_substring(qepcad,"(A v")==
        source.variables.size()-source.outer_variables.size(),
        "every pair-local source coordinate becomes one universal prenex variable");

    auto bad=p;
    bad.source_program_digest=std::string(64,'0');
    check(
        verify_qepcad_rcf_program(bad).status==
        QepcadRcfValidationStatus::InvalidSourceBinding,
        "backend verifier rejects Phase 6B2B source-digest tampering");

    bad=p;
    bad.maximum_total_degree+=1;
    check(
        verify_qepcad_rcf_program(bad).status==
        QepcadRcfValidationStatus::DegreeMismatch,
        "backend verifier rejects degree metadata tampering");

    if (!p.trig_definitions.empty()) {
        bad=p;
        bad.trig_definitions.front().cosine_lower=
            bad.trig_definitions.front().cosine_upper;
        check(
            verify_qepcad_rcf_program(bad).status==
            QepcadRcfValidationStatus::InvalidAlgebraicIsolation,
            "backend verifier rejects corrupted algebraic root isolation");
    }

    bool rejected_large_denominator=false;
    try {
        NetworkClosureInput large_den{
            2,{0,1},
            {
                Turn{PiRational{1,257}},
                Turn{PiRational{1,2}}
            }
        };
        (void)compile_qepcad_rcf_program(large_den);
    } catch (const std::invalid_argument&) {
        rejected_large_denominator=true;
    }
    check(
        rejected_large_denominator,
        "QEPCAD adapter fails closed beyond its explicit trig denominator resource limit");

    if (failures) {
        std::cerr << failures << " QEPCAD backend adapter test(s) failed\n";
        return 1;
    }

    std::cout << "QEPCAD exact RCF backend adapter tests passed\n";
    return 0;
}
