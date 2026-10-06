#pragma once

#include <pcg/semialgebraic_polynomial.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace pcg {

enum class QeVariableRole {
    AlgebraicCosine,
    AlgebraicSine,
    AlgebraicOperation,
    SourceOuter,
    SourcePairLocal
};

struct QeVariable {
    std::size_t index{};
    std::string name;
    QeVariableRole role{QeVariableRole::AlgebraicOperation};
    PiRational phase_pi;
    std::size_t source_variable{};
    std::size_t pair_e{};
    std::size_t pair_f{};
};

struct IntegerPolynomial {
    std::vector<BigInt> coefficients; // ascending degree order
};

struct TrigIsolationDefinition {
    PiRational phase_pi;
    std::size_t cosine_variable{};
    std::size_t sine_variable{};
    IntegerPolynomial cosine_polynomial;
    Rational cosine_lower{BigInt{-1}};
    Rational cosine_upper{BigInt{1}};
    int sine_sign{};
};

struct AlgebraicOperationDefinition {
    std::size_t variable{};
    std::string source_expression;
};

enum class RcfPolynomialKind {
    Rational,
    Variable,
    Add,
    Multiply,
    Negate
};

struct RcfPolynomialExpr {
    RcfPolynomialKind kind{RcfPolynomialKind::Rational};
    Rational rational{BigInt{0}};
    std::size_t variable{};
    std::vector<RcfPolynomialExpr> args;
};

enum class RcfRelation {
    EqualZero,
    GreaterEqualZero,
    GreaterZero
};

struct RcfAtom {
    RcfRelation relation{RcfRelation::EqualZero};
    RcfPolynomialExpr polynomial;
};

enum class RcfFormulaKind {
    True,
    False,
    Atom,
    And,
    Or,
    Not,
    Exists
};

struct RcfFormula {
    RcfFormulaKind kind{RcfFormulaKind::True};
    RcfAtom atom;
    std::vector<std::size_t> quantified_variables;
    std::vector<RcfFormula> children;
};

struct QepcadRcfProgram {
    std::string schema{"pcg-qepcad-rcf-program"};
    std::string schema_version{"1.0"};
    std::string backend_contract{"PCG-FT-QE-QEPCAD-001"};
    std::string theorem_id{"PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001"};
    SourceDigest source_digest{
        "SHA-256",
        "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208"
    };

    NetworkClosureInput input;
    std::string source_program_digest;

    std::vector<QeVariable> variables;
    std::vector<std::size_t> coefficient_variables;
    std::vector<TrigIsolationDefinition> trig_definitions;
    std::vector<AlgebraicOperationDefinition> operation_definitions;

    RcfFormula formula;

    std::size_t maximum_total_degree{};
    std::size_t rational_atom_count{};
};

enum class QepcadRcfValidationStatus {
    Valid,
    InvalidMetadata,
    InvalidSourceBinding,
    InvalidAlgebraicIsolation,
    InvalidStructure,
    DegreeMismatch
};

struct QepcadRcfValidation {
    QepcadRcfValidationStatus status{QepcadRcfValidationStatus::InvalidStructure};
    std::string detail;
};

QepcadRcfProgram lower_to_qepcad_rcf_program(
    const ExactSemialgebraicProgram& source);

QepcadRcfProgram compile_qepcad_rcf_program(
    const NetworkClosureInput& input);

std::size_t rcf_polynomial_degree(
    const RcfPolynomialExpr& expr);

std::string serialize_rcf_polynomial_expr(
    const RcfPolynomialExpr& expr);

std::string serialize_rcf_formula(
    const RcfFormula& formula);

std::string serialize_qepcad_rcf_program(
    const QepcadRcfProgram& program);

std::string qepcad_rcf_program_digest(
    const QepcadRcfProgram& program);

std::string export_qepcad_input(
    const QepcadRcfProgram& program);

QepcadRcfValidation verify_qepcad_rcf_program(
    const QepcadRcfProgram& program);

} // namespace pcg
