#include <pcg/network_closure.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace {

int failures=0;

void check(bool cond,const std::string& name) {
    if (!cond) {
        std::cerr << "FAIL: " << name << "\n";
        ++failures;
    }
}

pcg::NetworkClosureInput triangle(pcg::PiRational q) {
    return pcg::NetworkClosureInput{
        3,
        {0,1,2},
        {pcg::Turn{q},pcg::Turn{q},pcg::Turn{q}}
    };
}

} // namespace

int main() {
    using namespace pcg;

    // Three 120-degree turns produce chord directions at 60,180,300 degrees:
    // the unit positive kernel closes exactly.
    auto closed_tri=triangle(PiRational{2,3});
    auto c=solve_network_closure(closed_tri);
    check(c.status==NetworkClosureStatus::Closed,
          "equilateral direction triangle is exactly closed");
    check(c.proof==NetworkClosureProof::ExactPositiveRigidKernel,
          "closed triangle has exact positive rigid-kernel proof");
    check(c.assurance==ArithmeticAssurance::Exact,
          "closed triangle is exact");
    check(verify_network_closure_certificate(c.certificate).status==NetworkClosureStatus::Closed,
          "closed triangle certificate verifies");

    // Three 60-degree turns give directions 30,90,150 degrees, all in the
    // upper half plane. Their one-dimensional kernel cannot be strictly positive.
    auto open_tri=triangle(PiRational{1,3});
    auto o=solve_network_closure(open_tri);
    check(o.status==NetworkClosureStatus::NotClosed,
          "upper-half-plane triangle has no positive network closure");
    check(o.proof==NetworkClosureProof::ExactRigidSignObstruction,
          "nonclosed triangle has exact rigid sign obstruction");
    check(verify_network_closure_certificate(o.certificate).status==NetworkClosureStatus::NotClosed,
          "nonclosed triangle certificate verifies");

    // One quotient vertex and one self-loop imposes a nonzero directed chord
    // on a one-edge cycle; the exact matrix has full column rank.
    NetworkClosureInput loop{
        1,{0},{Turn{PiRational{1,2}}}
    };
    auto l=solve_network_closure(loop);
    check(l.status==NetworkClosureStatus::NotClosed,
          "single nonzero loop chord cannot close");
    check(l.proof==NetworkClosureProof::ExactFullRankObstruction,
          "single loop is rejected by exact full rank");

    // Production ABCABC rational witness already used by the fixed-turn solver.
    NetworkClosureInput abcabc{
        3,
        {0,1,2,0,1,2},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}}
        }
    };
    const auto basis=build_network_cycle_basis(abcabc);
    check(basis.rows.size()==4,
          "ABCABC quotient cycle-space dimension beta=4");
    auto a=solve_network_closure(abcabc);
    check(a.status==NetworkClosureStatus::Closed,
          "ABCABC production witness passes exact network closure");
    check(a.proof==NetworkClosureProof::ExactPositiveRigidKernel,
          "ABCABC witness is projectively rigid with positive kernel");
    check(verify_network_closure_certificate(a.certificate).status==NetworkClosureStatus::Closed,
          "ABCABC network certificate verifies");

    // Every cycle-basis row must have zero incidence at every quotient vertex.
    bool cycles_ok=true;
    for (const auto& row : basis.rows) {
        std::vector<int> div(3,0);
        for (std::size_t e=0;e<row.size();++e) {
            const std::size_t u=abcabc.trace_vertices[e];
            const std::size_t v=abcabc.trace_vertices[(e+1)%row.size()];
            div[u]+=row[e];
            div[v]-=row[e];
        }
        for (int x : div) if (x!=0) cycles_ok=false;
    }
    check(cycles_ok,"fundamental cycle basis lies in incidence kernel");

    // Cyclic reindexing changes the canonical input but must preserve closure.
    NetworkClosureInput rotated{
        3,
        {1,2,0,1,2,0},
        {
            Turn{PiRational{1,1}},Turn{PiRational{1,3}},
            Turn{PiRational{1,1}},Turn{PiRational{1,3}},
            Turn{PiRational{1,1}},Turn{PiRational{1,3}}
        }
    };
    check(solve_network_closure(rotated).status==NetworkClosureStatus::Closed,
          "cyclic source reindexing preserves network closure");

    // Software-domain cap: theorem remains general, but very high root-of-unity
    // degree is intentionally returned as INDETERMINATE rather than approximated.
    NetworkClosureInput high_order{
        2,{0,1},
        {Turn{PiRational{1,257}},Turn{PiRational{1,3}}}
    };
    auto h=solve_network_closure(high_order);
    check(h.status==NetworkClosureStatus::Indeterminate &&
          h.proof==NetworkClosureProof::CyclotomicOrderLimit,
          "high cyclotomic order is conservatively indeterminate");

    // Tamper resistance: exact positive-kernel coefficients are certificate evidence.
    auto tampered=a.certificate;
    check(!tampered.rigid_kernel.empty() &&
          !tampered.rigid_kernel[0].coefficients.empty(),
          "ABCABC certificate exposes exact algebraic kernel");
    tampered.rigid_kernel[0].coefficients[0]+=Rational(BigInt{1});
    auto bad=verify_network_closure_certificate(tampered);
    check(bad.status==NetworkClosureStatus::Indeterminate &&
          bad.termination==TerminationReason::BackendFailure,
          "verifier rejects tampered algebraic kernel");

    auto wrong_source=a.certificate;
    wrong_source.source_digest.value=std::string(64,'0');
    check(verify_network_closure_certificate(wrong_source).status==NetworkClosureStatus::Indeterminate,
          "verifier rejects tampered Network Closure source digest");

    auto wrong_rank=a.certificate;
    ++wrong_rank.exact_rank;
    check(verify_network_closure_certificate(wrong_rank).status==NetworkClosureStatus::Indeterminate,
          "verifier rejects tampered exact rank");

    if (failures) {
        std::cerr << failures << " network-closure test(s) failed\n";
        return 1;
    }

    std::cout << "Exact projectively-rigid Network Closure tests passed\n";
    return 0;
}
