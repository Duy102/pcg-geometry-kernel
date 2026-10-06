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

pcg::NetworkClosureInput n_cycle(std::size_t n, pcg::PiRational q) {
    pcg::NetworkClosureInput in;
    in.vertex_count=n;
    for (std::size_t i=0;i<n;++i) {
        in.trace_vertices.push_back(i);
        in.turns.emplace_back(q);
    }
    return in;
}

pcg::NetworkClosureInput four_cycle(pcg::PiRational q) {
    return n_cycle(4,q);
}

} // namespace

int main() {
    using namespace pcg;

    // Rigid closed triangle: unit positive kernel.
    auto closed_tri=triangle(PiRational{2,3});
    auto c=solve_network_closure(closed_tri);
    check(c.status==NetworkClosureStatus::Closed,
          "equilateral direction triangle is closed");
    check(c.proof==NetworkClosureProof::CertifiedPositiveKernel,
          "closed triangle has certified positive-kernel proof");
    check(c.assurance==ArithmeticAssurance::CertifiedNumerical,
          "positive-kernel conclusion uses exact algebra plus certified sign separation");
    check(c.certificate.schema_version=="1.1",
          "general Network Closure certificate schema is v1.1");
    check(verify_network_closure_certificate(c.certificate).status==NetworkClosureStatus::Closed,
          "closed triangle certificate verifies");

    // Rigid non-closed triangle: all directions lie in the upper half plane.
    auto open_tri=triangle(PiRational{1,3});
    auto o=solve_network_closure(open_tri);
    check(o.status==NetworkClosureStatus::NotClosed,
          "upper-half-plane triangle has no positive network closure");
    check(o.proof==NetworkClosureProof::CertifiedStiemkeObstruction,
          "nonclosed triangle has certified Stiemke obstruction");
    check(o.assurance==ArithmeticAssurance::CertifiedNumerical,
          "Stiemke obstruction uses exact algebra plus certified sign separation");
    check(!o.certificate.stiemke_dual.empty(),
          "Stiemke obstruction stores a dual witness");
    check(verify_network_closure_certificate(o.certificate).status==NetworkClosureStatus::NotClosed,
          "nonclosed triangle certificate verifies");

    // One quotient vertex and one self-loop gives full column rank.
    NetworkClosureInput loop{
        1,{0},{Turn{PiRational{1,2}}}
    };
    auto l=solve_network_closure(loop);
    check(l.status==NetworkClosureStatus::NotClosed,
          "single nonzero loop chord cannot close");
    check(l.proof==NetworkClosureProof::ExactFullRankObstruction,
          "single loop is rejected by exact full rank");
    check(l.assurance==ArithmeticAssurance::Exact,
          "full-rank obstruction is exact without sign numerics");

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
          "ABCABC production witness passes Network Closure");
    check(a.proof==NetworkClosureProof::CertifiedPositiveKernel,
          "ABCABC witness has a certified positive kernel");
    check(verify_network_closure_certificate(a.certificate).status==NetworkClosureStatus::Closed,
          "ABCABC network certificate verifies");

    // Every fundamental cycle row must lie in the incidence kernel.
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

    // Phase 5B target #1: beta=1, rank=2, nullity=2, with positive kernel.
    // Four half-turns give chord directions 45,135,225,315 degrees.
    auto square_closed=four_cycle(PiRational{1,2});
    auto sc=solve_network_closure(square_closed);
    check(sc.certificate.exact_rank==2,
          "closed four-cycle has exact rank two");
    check(square_closed.turns.size()-sc.certificate.exact_rank==2,
          "closed four-cycle has higher-dimensional kernel");
    check(sc.status==NetworkClosureStatus::Closed,
          "higher-dimensional square kernel contains a positive vector");
    check(sc.proof==NetworkClosureProof::CertifiedPositiveKernel,
          "higher-dimensional closure uses general positive feasibility");
    check(sc.certificate.positive_kernel.size()==4,
          "higher-dimensional closure carries a four-component positive witness");
    check(verify_network_closure_certificate(sc.certificate).status==NetworkClosureStatus::Closed,
          "higher-dimensional positive-kernel certificate verifies");

    // Phase 5B target #2: same beta/rank/nullity, but every direction has
    // positive y component, so no strictly positive chord combination can sum to zero.
    // Stiemke supplies a nonzero dual stress N^T y >= 0.
    auto square_open=four_cycle(PiRational{1,4});
    auto so=solve_network_closure(square_open);
    check(so.certificate.exact_rank==2,
          "open four-cycle has exact rank two");
    check(square_open.turns.size()-so.certificate.exact_rank==2,
          "open four-cycle also has higher-dimensional kernel");
    check(so.status==NetworkClosureStatus::NotClosed,
          "higher-dimensional upper-half-plane directions are not closed");
    check(so.proof==NetworkClosureProof::CertifiedStiemkeObstruction,
          "higher-dimensional nonclosure carries a Stiemke dual");
    check(!so.certificate.stiemke_dual.empty(),
          "higher-dimensional Stiemke dual is stored");
    check(verify_network_closure_certificate(so.certificate).status==NetworkClosureStatus::NotClosed,
          "higher-dimensional Stiemke certificate verifies");

    // General higher-dimensional family:
    // tau=2*pi/n gives n evenly-spaced chord directions around the full circle,
    // so c=(1,...,1) is positive. tau=pi/n keeps all chord directions in one
    // open half-plane, so Stiemke must certify non-closure.
    for (std::size_t n=4;n<=8;++n) {
        auto regular=n_cycle(n,PiRational{BigInt{2},BigInt{n}});
        auto rc=solve_network_closure(regular);
        check(rc.certificate.exact_rank==2,
              "regular n-cycle has exact rank two");
        check(n-rc.certificate.exact_rank==n-2,
              "regular n-cycle exercises expected higher nullity");
        check(rc.status==NetworkClosureStatus::Closed &&
              rc.proof==NetworkClosureProof::CertifiedPositiveKernel,
              "regular n-cycle has certified positive kernel");
        check(verify_network_closure_certificate(rc.certificate).status==
              NetworkClosureStatus::Closed,
              "regular n-cycle certificate verifies");

        auto halfplane=n_cycle(n,PiRational{BigInt{1},BigInt{n}});
        auto hp=solve_network_closure(halfplane);
        check(hp.certificate.exact_rank==2,
              "half-plane n-cycle has exact rank two");
        check(n-hp.certificate.exact_rank==n-2,
              "half-plane n-cycle exercises expected higher nullity");
        check(hp.status==NetworkClosureStatus::NotClosed &&
              hp.proof==NetworkClosureProof::CertifiedStiemkeObstruction,
              "half-plane n-cycle has certified Stiemke obstruction");
        check(verify_network_closure_certificate(hp.certificate).status==
              NetworkClosureStatus::NotClosed,
              "half-plane n-cycle Stiemke certificate verifies");
    }

    // Cyclic source reindexing preserves the geometric closure problem.
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
    // order is returned conservatively as INDETERMINATE.
    NetworkClosureInput high_order{
        2,{0,1},
        {Turn{PiRational{1,257}},Turn{PiRational{1,3}}}
    };
    auto h=solve_network_closure(high_order);
    check(h.status==NetworkClosureStatus::Indeterminate &&
          h.proof==NetworkClosureProof::CyclotomicOrderLimit,
          "high cyclotomic order is conservatively indeterminate");

    // Tamper resistance for positive-kernel evidence.
    auto tampered=a.certificate;
    check(!tampered.positive_kernel.empty() &&
          !tampered.positive_kernel[0].coefficients.empty(),
          "ABCABC certificate exposes algebraic positive kernel");
    tampered.positive_kernel[0].coefficients[0]+=Rational(BigInt{1});
    auto bad=verify_network_closure_certificate(tampered);
    check(bad.status==NetworkClosureStatus::Indeterminate &&
          bad.termination==TerminationReason::BackendFailure,
          "verifier rejects tampered positive kernel");

    // Tamper resistance for Stiemke dual evidence.
    auto tampered_dual=so.certificate;
    check(!tampered_dual.stiemke_dual.empty() &&
          !tampered_dual.stiemke_dual[0].coefficients.empty(),
          "nonclosed certificate exposes algebraic Stiemke dual");
    tampered_dual.stiemke_dual[0].coefficients[0]+=Rational(BigInt{1});
    auto bad_dual=verify_network_closure_certificate(tampered_dual);
    check(bad_dual.status==NetworkClosureStatus::Indeterminate &&
          bad_dual.termination==TerminationReason::BackendFailure,
          "verifier rejects tampered Stiemke dual");

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

    std::cout << "General certified Network Closure tests passed\n";
    return 0;
}
