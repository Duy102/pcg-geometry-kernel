#pragma once

#include <pcg/semialgebraic.hpp>

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

namespace pcg {

enum class ExactAlgebraicKind {
    Rational,
    SinPi,
    CosPi,
    Add,
    Multiply,
    Negate,
    Inverse
};

struct ExactAlgebraicExpr {
    ExactAlgebraicKind kind{ExactAlgebraicKind::Rational};
    Rational rational{BigInt{0}};
    PiRational phase_pi;
    std::vector<ExactAlgebraicExpr> args;
};

enum class PolynomialExprKind {
    Constant,
    Variable,
    Add,
    Multiply,
    Negate
};

struct PolynomialExpr {
    PolynomialExprKind kind{PolynomialExprKind::Constant};
    ExactAlgebraicExpr constant;
    std::size_t variable{};
    std::vector<PolynomialExpr> args;
};

enum class PolynomialRelation {
    EqualZero,
    GreaterEqualZero,
    GreaterZero
};

struct PolynomialAtom {
    PolynomialRelation relation{PolynomialRelation::EqualZero};
    PolynomialExpr polynomial;
};

enum class FormulaKind {
    True,
    False,
    Atom,
    And,
    Or,
    Not,
    Exists
};

struct SemialgebraicFormula {
    FormulaKind kind{FormulaKind::True};
    PolynomialAtom atom;
    std::vector<std::size_t> quantified_variables;
    std::vector<SemialgebraicFormula> children;
};

enum class SemialgebraicVariableScope {
    Outer,
    PairLocal
};

struct SemialgebraicVariable {
    std::size_t index{};
    std::string name;
    SemialgebraicVariableScope scope{SemialgebraicVariableScope::Outer};
    std::size_t pair_e{};
    std::size_t pair_f{};
};

struct ExactSemialgebraicProgram {
    std::string schema{"pcg-exact-semialgebraic-program"};
    std::string schema_version{"1.0"};
    std::string lowering_contract{"PCG-FT-SA-POLY-AST-001"};
    std::string theorem_id{"PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001"};
    SourceDigest source_digest{
        "SHA-256",
        "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208"
    };

    NetworkClosureInput input;
    std::string source_ir_digest;

    std::vector<SemialgebraicVariable> variables;
    std::vector<std::size_t> outer_variables;
    SemialgebraicFormula formula;

    std::size_t maximum_polynomial_degree{};
    std::size_t polynomial_atom_count{};
};

enum class ExactSemialgebraicProgramValidationStatus {
    Valid,
    InvalidMetadata,
    InvalidSourceIrBinding,
    InvalidStructure,
    DegreeContractViolation
};

struct ExactSemialgebraicProgramValidation {
    ExactSemialgebraicProgramValidationStatus status{
        ExactSemialgebraicProgramValidationStatus::InvalidStructure};
    std::string detail;
};

ExactSemialgebraicProgram lower_fixed_turn_semialgebraic_problem(
    const FixedTurnSemialgebraicProblem& problem);

ExactSemialgebraicProgram compile_exact_semialgebraic_program(
    const NetworkClosureInput& input);

std::string serialize_exact_algebraic_expr(
    const ExactAlgebraicExpr& expr);

std::string serialize_polynomial_expr(
    const PolynomialExpr& expr);

std::string serialize_semialgebraic_formula(
    const SemialgebraicFormula& formula);

std::string serialize_exact_semialgebraic_program(
    const ExactSemialgebraicProgram& program);

std::string exact_semialgebraic_program_digest(
    const ExactSemialgebraicProgram& program);

std::size_t polynomial_degree(
    const PolynomialExpr& expr);

ExactSemialgebraicProgramValidation verify_exact_semialgebraic_program(
    const ExactSemialgebraicProgram& program);

} // namespace pcg
