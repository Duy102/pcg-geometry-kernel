#include <pcg/post_closure.hpp>

#include <iostream>
#include <string>

namespace {
int failures=0;
void check(bool cond,const std::string& name) {
    if (!cond) {
        std::cerr << "FAIL: " << name << "\n";
        ++failures;
    }
}
}

int main() {
    using namespace pcg;

    // General remote Hit: two distinct unit semicircles intersect in the
    // interior of both finite arcs.
    CertifiedFiniteArc upper{
        {CertifiedScalar(1.0,1.0),CertifiedScalar(0.0,0.0)},
        {CertifiedScalar(-1.0,-1.0),CertifiedScalar(0.0,0.0)},
        PiRational{1,2},Turn{PiRational{1,1}},PositiveInterval::point(2.0)
    };
    CertifiedFiniteArc lower_shifted{
        {CertifiedScalar(-1.0,-1.0),CertifiedScalar(0.5,0.5)},
        {CertifiedScalar(1.0,1.0),CertifiedScalar(0.5,0.5)},
        PiRational{-1,2},Turn{PiRational{1,1}},PositiveInterval::point(2.0)
    };
    check(pcg_hit(upper,lower_shifted)==CertifiedTruth::True,
          "general finite-arc Hit certifies an interior remote intersection");

    CertifiedFiniteArc far_lower{
        {CertifiedScalar(-1.0,-1.0),CertifiedScalar(3.0,3.0)},
        {CertifiedScalar(1.0,1.0),CertifiedScalar(3.0,3.0)},
        PiRational{-1,2},Turn{PiRational{1,1}},PositiveInterval::point(2.0)
    };
    check(pcg_hit(upper,far_lower)==CertifiedTruth::False,
          "general finite-arc Hit certifies separated supporting circles");

    // C1-adjacent arcs share an allowed endpoint with one tangent line.
    // Distinct tangent support circles have no second intersection.
    SharedStartArc tangent_a{
        PiRational{0,1},Turn{PiRational{1,3}},PositiveInterval::point(1.0)};
    SharedStartArc tangent_b{
        PiRational{1,1},Turn{PiRational{-1,2}},PositiveInterval::point(1.0)};
    check(sin_pi_sign(tangent_b.tangent_phase_pi-tangent_a.tangent_phase_pi)==0,
          "tangent-safe regression uses one exact tangent line");
    check(shared_start_support_relation(tangent_a,tangent_b)==SupportRelation::Distinct,
          "tangent-safe regression has distinct support circles");

    // Production ABCABC witness: Phase 5B closure is projectively rigid and
    // Phase 6A must recover the same trace-faithful conclusion as the
    // specialized theorem path.
    NetworkClosureInput abcabc{
        3,{0,1,2,0,1,2},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}}
        }
    };
    auto nc=solve_network_closure(abcabc);
    check(nc.status==NetworkClosureStatus::Closed,
          "ABCABC witness closes before post-closure geometry");
    check(abcabc.turns.size()-nc.certificate.exact_rank==1,
          "ABCABC witness is projectively rigid");

    auto emb=reconstruct_projectively_rigid_embedding(nc.certificate);
    check(emb.status==ProjectivelyRigidEmbeddingStatus::Ready,
          "projectively rigid embedding reconstructs exactly");
    check(emb.embedding.vertices.size()==3 &&
          emb.embedding.chord_magnitudes.size()==6,
          "embedding exports quotient vertices and normalized chord magnitudes");

    auto post=solve_projectively_rigid_post_closure(nc.certificate);
    check(post.decision==Decision::Realizable,
          "ABCABC witness passes general projectively-rigid post-closure verifier");
    check(post.proof==RigidPostClosureProof::AllPairsClear,
          "trace-faithful witness has all-pairs-clear proof");
    check(post.certificate.checked_pairs==15,
          "six finite arcs check all fifteen unordered pairs");
    check(verify_projectively_rigid_post_closure_certificate(post.certificate).decision==
          Decision::Realizable,
          "post-closure certificate verifies");

    // Phase 6B1: a higher-dimensional positive-kernel witness may itself be
    // enough to certify existence, even though it is not a complete search of
    // the positive closure polytope.
    NetworkClosureInput convex_quad{
        4,{0,1,2,3},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,2}},
            Turn{PiRational{1,4}},Turn{PiRational{11,12}}
        }
    };
    auto cq=solve_network_closure(convex_quad);
    check(cq.status==NetworkClosureStatus::Closed &&
          convex_quad.turns.size()-cq.certificate.exact_rank==2,
          "generic convex quadrilateral has two-dimensional positive closure");
    auto cqemb=reconstruct_positive_kernel_embedding(cq.certificate);
    check(cqemb.status==ProjectivelyRigidEmbeddingStatus::Ready,
          "higher-dimensional positive-kernel witness reconstructs exactly");
    auto cqpost=solve_positive_kernel_post_closure_witness(cq.certificate);
    check(cqpost.decision==Decision::Realizable &&
          cqpost.proof==PositiveKernelPostClosureProof::TraceFaithfulMetricWitness,
          "higher-dimensional positive-kernel witness certifies realizability");
    check(cqpost.certificate.checked_pairs==6,
          "four finite arcs check all six unordered pairs");
    check(verify_positive_kernel_post_closure_certificate(cqpost.certificate).decision==
          Decision::Realizable,
          "higher-dimensional witness certificate verifies");

    auto cqbad=cqpost.certificate;
    cqbad.source_digest.value=std::string(64,'0');
    check(verify_positive_kernel_post_closure_certificate(cqbad).termination==
          TerminationReason::BackendFailure,
          "higher-dimensional witness verifier rejects theorem-source tampering");

    // Phase 6A still deliberately rejects higher-dimensional closure because
    // its NOT_REALIZABLE conclusions rely on uniqueness of the normalized
    // metric skeleton.
    NetworkClosureInput square{
        4,{0,1,2,3},
        {
            Turn{PiRational{1,2}},Turn{PiRational{1,2}},
            Turn{PiRational{1,2}},Turn{PiRational{1,2}}
        }
    };
    auto sq=solve_network_closure(square);
    check(sq.status==NetworkClosureStatus::Closed &&
          square.turns.size()-sq.certificate.exact_rank==2,
          "square exercises higher-dimensional positive closure");
    auto sqpost=solve_projectively_rigid_post_closure(sq.certificate);
    check(sqpost.decision==Decision::Indeterminate &&
          sqpost.proof==RigidPostClosureProof::HigherDimensionalClosure,
          "Phase 6A does not overclaim higher-dimensional length selection");

    auto sqwitness=solve_positive_kernel_post_closure_witness(sq.certificate);
    check(sqwitness.decision!=Decision::NotRealizable,
          "one higher-dimensional metric witness can never certify global non-realizability");

    // Upstream nonclosure propagates as a realizability obstruction.
    NetworkClosureInput open_tri{
        3,{0,1,2},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,3}},Turn{PiRational{1,3}}
        }
    };
    auto open=solve_network_closure(open_tri);
    auto open_post=solve_projectively_rigid_post_closure(open.certificate);
    check(open_post.decision==Decision::NotRealizable &&
          open_post.proof==RigidPostClosureProof::NetworkClosureObstruction,
          "Network Closure obstruction propagates through Phase 6A");
    auto open_witness=solve_positive_kernel_post_closure_witness(open.certificate);
    check(open_witness.decision==Decision::NotRealizable &&
          open_witness.proof==PositiveKernelPostClosureProof::NetworkClosureObstruction,
          "global Network Closure obstruction propagates through Phase 6B1");

    auto tampered=post.certificate;
    tampered.source_digest.value=std::string(64,'0');
    auto bad=verify_projectively_rigid_post_closure_certificate(tampered);
    check(bad.decision==Decision::Indeterminate &&
          bad.termination==TerminationReason::BackendFailure,
          "post-closure verifier rejects theorem-source tampering");

    if (failures) {
        std::cerr << failures << " post-closure test(s) failed\n";
        return 1;
    }
    std::cout << "Post-closure geometry tests passed\n";
    return 0;
}
