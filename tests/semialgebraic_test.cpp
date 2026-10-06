#include <pcg/semialgebraic.hpp>

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

    NetworkClosureInput mixed_triangle{
        3,{0,1,2},
        {
            Turn{PiRational{1,2}},
            Turn{PiRational{1,1}},
            Turn{PiRational{3,2}}
        }
    };
    auto p=compile_fixed_turn_semialgebraic_problem(mixed_triangle);

    check(p.edges.size()==3,
          "compiler emits one theorem edge specification per source edge");
    check(p.edges[0].arc_class==SemialgebraicArcClass::Minor &&
          p.edges[1].arc_class==SemialgebraicArcClass::Semicircle &&
          p.edges[2].arc_class==SemialgebraicArcClass::Major,
          "compiler preserves minor, semicircle, and major membership branches");
    check(p.distinct_vertex_pairs.size()==3,
          "three quotient vertices emit all three distinctness constraints");
    check(p.forbidden_pairs.size()==3,
          "three finite arcs emit all three unordered bad-pair formulas");
    check(p.forbidden_pairs[0].allowed_common_vertices.size()==1,
          "adjacent triangle arcs retain their allowed common quotient vertex");
    check(p.outer_real_variable_count==9 &&
          p.pair_local_real_variable_count==6,
          "quantifier accounting matches 2V+m outer and two q coordinates per pair");
    check(verify_fixed_turn_semialgebraic_problem(p).status==
          SemialgebraicProblemValidationStatus::Valid,
          "deterministic theorem IR verifies");

    const auto d1=fixed_turn_semialgebraic_problem_digest(p);
    const auto d2=fixed_turn_semialgebraic_problem_digest(
        compile_fixed_turn_semialgebraic_problem(mixed_triangle));
    check(d1==d2 && d1.size()==64,
          "canonical semialgebraic problem digest is deterministic");

    auto tampered=p;
    tampered.forbidden_pairs.pop_back();
    check(verify_fixed_turn_semialgebraic_problem(tampered).status==
          SemialgebraicProblemValidationStatus::InvalidStructure,
          "verifier rejects missing forbidden pair");

    tampered=p;
    tampered.source_digest.value=std::string(64,'0');
    check(verify_fixed_turn_semialgebraic_problem(tampered).status==
          SemialgebraicProblemValidationStatus::InvalidMetadata,
          "verifier rejects theorem-source tampering");

    NetworkClosureInput abcabc{
        3,{0,1,2,0,1,2},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}}
        }
    };
    auto a=compile_fixed_turn_semialgebraic_problem(abcabc);
    check(a.forbidden_pairs.size()==15,
          "six arcs compile all fifteen pairwise incidence formulas");
    check(a.prescribed_transversality.size()==3,
          "ABCABC compiles one fixed tangent-transversality check per repeated symbol");
    bool all_transverse=true;
    for (const auto& x:a.prescribed_transversality)
        all_transverse=all_transverse && x.transverse;
    check(all_transverse,
          "ABCABC prescribed double points are transverse in fixed turn data");

    NetworkClosureInput convex_quad{
        4,{0,1,2,3},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,2}},
            Turn{PiRational{1,4}},Turn{PiRational{11,12}}
        }
    };
    auto q=compile_fixed_turn_semialgebraic_problem(convex_quad);
    check(q.distinct_vertex_pairs.size()==6 &&
          q.forbidden_pairs.size()==6,
          "Phase 6B1 nullity-two example compiles the full finite formula skeleton");

    if (failures) {
        std::cerr << failures << " semialgebraic IR test(s) failed\n";
        return 1;
    }
    std::cout << "Semialgebraic IR tests passed\n";
    return 0;
}
