#include <pcg/semialgebraic_polynomial.hpp>

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

std::size_t count_kind(
    const pcg::SemialgebraicFormula& f,
    pcg::FormulaKind kind) {
    std::size_t n=f.kind==kind?1:0;
    for (const auto& c:f.children) n+=count_kind(c,kind);
    return n;
}

bool all_atoms_quadratic_or_less(
    const pcg::SemialgebraicFormula& f) {
    if (f.kind==pcg::FormulaKind::Atom &&
        pcg::polynomial_degree(f.atom.polynomial)>2)
        return false;
    for (const auto& c:f.children)
        if (!all_atoms_quadratic_or_less(c)) return false;
    return true;
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

    const auto ir=compile_fixed_turn_semialgebraic_problem(mixed_triangle);
    auto p=lower_fixed_turn_semialgebraic_problem(ir);

    check(p.source_ir_digest==
          fixed_turn_semialgebraic_problem_digest(ir),
          "Phase 6B2B binds exactly to the verified Phase 6B2A IR");
    check(p.variables.size()==15 && p.outer_variables.size()==9,
          "triangle lowers 9 outer variables plus six pair-local q coordinates");
    check(p.formula.kind==FormulaKind::Exists &&
          p.formula.quantified_variables.size()==9,
          "complete theorem sentence starts with existential p,c quantification");
    check(count_kind(p.formula,FormulaKind::Exists)==4,
          "triangle has one outer exists plus one local q exists for each of three pairs");
    check(count_kind(p.formula,FormulaKind::Not)==3,
          "every unordered arc pair lowers to not exists q I_ef");
    check(p.maximum_polynomial_degree==2 &&
          all_atoms_quadratic_or_less(p.formula),
          "all lowered geometric atoms stay within the theorem's quadratic degree bound");
    check(p.polynomial_atom_count>0,
          "exact polynomial AST contains explicit atomic constraints");

    const auto serialized=serialize_exact_semialgebraic_program(p);
    check(serialized.find("sin_pi")!=std::string::npos &&
          serialized.find("cos_pi")!=std::string::npos &&
          serialized.find("\"inv\"")!=std::string::npos,
          "coefficients preserve exact rational-pi algebraic trigonometric structure");
    check(serialized.find("0.333") == std::string::npos,
          "lowering introduces no decimal approximation for theorem coefficients");

    const auto d1=exact_semialgebraic_program_digest(p);
    const auto d2=exact_semialgebraic_program_digest(
        compile_exact_semialgebraic_program(mixed_triangle));
    check(d1==d2 && d1.size()==64,
          "exact polynomial program SHA-256 is deterministic");

    check(verify_exact_semialgebraic_program(p).status==
          ExactSemialgebraicProgramValidationStatus::Valid,
          "exact semialgebraic program recompiles and verifies");

    auto bad=p;
    bad.source_ir_digest=std::string(64,'0');
    check(verify_exact_semialgebraic_program(bad).status==
          ExactSemialgebraicProgramValidationStatus::InvalidSourceIrBinding,
          "verifier rejects Phase 6B2A binding tampering");

    bad=p;
    bad.maximum_polynomial_degree=3;
    check(verify_exact_semialgebraic_program(bad).status==
          ExactSemialgebraicProgramValidationStatus::DegreeContractViolation,
          "verifier rejects a mutated quadratic-degree contract");

    NetworkClosureInput abcabc{
        3,{0,1,2,0,1,2},
        {
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}},
            Turn{PiRational{1,3}},Turn{PiRational{1,1}}
        }
    };
    auto a=compile_exact_semialgebraic_program(abcabc);
    check(a.outer_variables.size()==12,
          "ABCABC has six vertex coordinates plus six positive chord variables");
    check(a.variables.size()==42,
          "ABCABC adds thirty local coordinates for its fifteen unordered arc pairs");
    check(count_kind(a.formula,FormulaKind::Not)==15 &&
          count_kind(a.formula,FormulaKind::Exists)==16,
          "ABCABC lowers all fifteen forbidden-pair existential incidences");

    // Fixed non-transversality is rejected at the formula level before metric
    // length search: repeated symbol 0 has tangent phases separated by pi.
    NetworkClosureInput nontransverse{
        2,{0,1,0,1},
        {
            Turn{PiRational{1,2}},Turn{PiRational{1,2}},
            Turn{PiRational{1,2}},Turn{PiRational{1,2}}
        }
    };
    auto nt=compile_exact_semialgebraic_program(nontransverse);
    check(count_kind(nt.formula,FormulaKind::False)>0,
          "fixed prescribed non-transversality makes the theorem formula false");

    if (failures) {
        std::cerr << failures << " exact semialgebraic lowering test(s) failed\n";
        return 1;
    }
    std::cout << "Exact semialgebraic polynomial AST tests passed\n";
    return 0;
}
