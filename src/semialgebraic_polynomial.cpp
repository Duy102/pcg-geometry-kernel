#include <pcg/semialgebraic_polynomial.hpp>

#include <algorithm>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pcg {
namespace {

constexpr const char* kSchema="pcg-exact-semialgebraic-program";
constexpr const char* kSchemaVersion="1.0";
constexpr const char* kLoweringContract="PCG-FT-SA-POLY-AST-001";
constexpr const char* kTheoremId="PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001";
constexpr const char* kSourceDigest=
    "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208";

std::string rational_text(const Rational& r) {
    std::ostringstream os;
    os << r.numerator() << "/" << r.denominator();
    return os.str();
}

ExactAlgebraicExpr a_rat(Rational q) {
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::Rational;
    x.rational=std::move(q);
    return x;
}

ExactAlgebraicExpr a_int(long long n) {
    return a_rat(Rational(BigInt{n}));
}

ExactAlgebraicExpr a_sin(PiRational p) {
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::SinPi;
    x.phase_pi=std::move(p);
    return x;
}

ExactAlgebraicExpr a_cos(PiRational p) {
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::CosPi;
    x.phase_pi=std::move(p);
    return x;
}

ExactAlgebraicExpr a_add(std::vector<ExactAlgebraicExpr> xs) {
    if (xs.empty()) return a_int(0);
    if (xs.size()==1) return std::move(xs.front());
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::Add;
    x.args=std::move(xs);
    return x;
}

ExactAlgebraicExpr a_mul(std::vector<ExactAlgebraicExpr> xs) {
    if (xs.empty()) return a_int(1);
    if (xs.size()==1) return std::move(xs.front());
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::Multiply;
    x.args=std::move(xs);
    return x;
}

ExactAlgebraicExpr a_neg(ExactAlgebraicExpr a) {
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::Negate;
    x.args.push_back(std::move(a));
    return x;
}

ExactAlgebraicExpr a_inv(ExactAlgebraicExpr a) {
    ExactAlgebraicExpr x;
    x.kind=ExactAlgebraicKind::Inverse;
    x.args.push_back(std::move(a));
    return x;
}

PolynomialExpr p_const(ExactAlgebraicExpr a) {
    PolynomialExpr p;
    p.kind=PolynomialExprKind::Constant;
    p.constant=std::move(a);
    return p;
}

PolynomialExpr p_int(long long n) {
    return p_const(a_int(n));
}

PolynomialExpr p_var(std::size_t v) {
    PolynomialExpr p;
    p.kind=PolynomialExprKind::Variable;
    p.variable=v;
    return p;
}

PolynomialExpr p_add(std::vector<PolynomialExpr> xs) {
    if (xs.empty()) return p_int(0);
    if (xs.size()==1) return std::move(xs.front());
    PolynomialExpr p;
    p.kind=PolynomialExprKind::Add;
    p.args=std::move(xs);
    return p;
}

PolynomialExpr p_mul(std::vector<PolynomialExpr> xs) {
    if (xs.empty()) return p_int(1);
    if (xs.size()==1) return std::move(xs.front());
    PolynomialExpr p;
    p.kind=PolynomialExprKind::Multiply;
    p.args=std::move(xs);
    return p;
}

PolynomialExpr p_neg(PolynomialExpr a) {
    PolynomialExpr p;
    p.kind=PolynomialExprKind::Negate;
    p.args.push_back(std::move(a));
    return p;
}

PolynomialExpr p_sub(PolynomialExpr a,PolynomialExpr b) {
    std::vector<PolynomialExpr> xs;
    xs.push_back(std::move(a));
    xs.push_back(p_neg(std::move(b)));
    return p_add(std::move(xs));
}

PolynomialExpr p_square(PolynomialExpr a) {
    PolynomialExpr b=a;
    std::vector<PolynomialExpr> xs;
    xs.push_back(std::move(a));
    xs.push_back(std::move(b));
    return p_mul(std::move(xs));
}

PolynomialExpr p_scale(ExactAlgebraicExpr a,PolynomialExpr p) {
    std::vector<PolynomialExpr> xs;
    xs.push_back(p_const(std::move(a)));
    xs.push_back(std::move(p));
    return p_mul(std::move(xs));
}

SemialgebraicFormula f_true() {
    SemialgebraicFormula f;
    f.kind=FormulaKind::True;
    return f;
}

SemialgebraicFormula f_false() {
    SemialgebraicFormula f;
    f.kind=FormulaKind::False;
    return f;
}

SemialgebraicFormula f_atom(PolynomialRelation rel,PolynomialExpr p) {
    SemialgebraicFormula f;
    f.kind=FormulaKind::Atom;
    f.atom.relation=rel;
    f.atom.polynomial=std::move(p);
    return f;
}

SemialgebraicFormula f_and(std::vector<SemialgebraicFormula> xs) {
    if (xs.empty()) return f_true();
    if (xs.size()==1) return std::move(xs.front());
    SemialgebraicFormula f;
    f.kind=FormulaKind::And;
    f.children=std::move(xs);
    return f;
}

SemialgebraicFormula f_or(std::vector<SemialgebraicFormula> xs) {
    if (xs.empty()) return f_false();
    if (xs.size()==1) return std::move(xs.front());
    SemialgebraicFormula f;
    f.kind=FormulaKind::Or;
    f.children=std::move(xs);
    return f;
}

SemialgebraicFormula f_not(SemialgebraicFormula x) {
    SemialgebraicFormula f;
    f.kind=FormulaKind::Not;
    f.children.push_back(std::move(x));
    return f;
}

SemialgebraicFormula f_exists(
    std::vector<std::size_t> vars,
    SemialgebraicFormula body) {
    SemialgebraicFormula f;
    f.kind=FormulaKind::Exists;
    f.quantified_variables=std::move(vars);
    f.children.push_back(std::move(body));
    return f;
}

struct EdgePolynomials {
    PolynomialExpr ox;
    PolynomialExpr oy;
    PolynomialExpr rminus_x;
    PolynomialExpr rminus_y;
    PolynomialExpr rplus_x;
    PolynomialExpr rplus_y;
    ExactAlgebraicExpr lambda;
};

std::size_t px(std::size_t v) { return 2*v; }
std::size_t py(std::size_t v) { return 2*v+1; }

std::size_t cvar(std::size_t vertex_count,std::size_t e) {
    return 2*vertex_count+e;
}

EdgePolynomials derive_edge(
    const SemialgebraicEdgeSpec& e,
    std::size_t vertex_count) {

    const auto ct=a_cos(e.tangent_phase_pi);
    const auto st=a_sin(e.tangent_phase_pi);
    const auto half_turn=e.turn_pi/2;
    const auto lambda=a_inv(a_mul({a_int(2),a_sin(half_turn)}));

    const auto c=p_var(cvar(vertex_count,e.edge));
    const auto ptx=p_var(px(e.tail));
    const auto pty=p_var(py(e.tail));
    const auto phx=p_var(px(e.head));
    const auto phy=p_var(py(e.head));

    const auto lambda_st=a_mul({lambda,st});
    const auto lambda_ct=a_mul({lambda,ct});

    // o = p_tail + lambda*c*J(t), J(tx,ty)=(-ty,tx).
    auto ox=p_sub(ptx,p_scale(lambda_st,c));
    auto oy=p_add({pty,p_scale(lambda_ct,c)});

    auto rminus_x=p_scale(lambda_st,c);
    auto rminus_y=p_neg(p_scale(lambda_ct,c));
    auto rplus_x=p_sub(phx,ox);
    auto rplus_y=p_sub(phy,oy);

    return {
        std::move(ox),std::move(oy),
        std::move(rminus_x),std::move(rminus_y),
        std::move(rplus_x),std::move(rplus_y),
        lambda
    };
}

SemialgebraicFormula arc_membership(
    const SemialgebraicEdgeSpec& e,
    const EdgePolynomials& d,
    std::size_t qx,
    std::size_t qy,
    std::size_t vertex_count) {

    const auto c=p_var(cvar(vertex_count,e.edge));
    auto x=p_sub(p_var(qx),d.ox);
    auto y=p_sub(p_var(qy),d.oy);

    auto circle=p_sub(
        p_add({p_square(x),p_square(y)}),
        p_square(p_scale(d.lambda,c)));

    auto cross_a=p_sub(
        p_mul({d.rminus_x,y}),
        p_mul({d.rminus_y,x}));
    if (e.turn_sign<0) cross_a=p_neg(std::move(cross_a));

    auto cross_b=p_sub(
        p_mul({x,d.rplus_y}),
        p_mul({y,d.rplus_x}));
    if (e.turn_sign<0) cross_b=p_neg(std::move(cross_b));

    auto f0=f_atom(PolynomialRelation::EqualZero,std::move(circle));
    auto age=f_atom(PolynomialRelation::GreaterEqualZero,std::move(cross_a));
    auto bge=f_atom(PolynomialRelation::GreaterEqualZero,std::move(cross_b));

    if (e.arc_class==SemialgebraicArcClass::Minor)
        return f_and({std::move(f0),std::move(age),std::move(bge)});
    if (e.arc_class==SemialgebraicArcClass::Semicircle)
        return f_and({std::move(f0),std::move(age)});
    return f_and({
        std::move(f0),
        f_or({std::move(age),std::move(bge)})
    });
}

SemialgebraicFormula distinct_point(
    std::size_t ax,std::size_t ay,
    std::size_t bx,std::size_t by) {
    auto dx=p_sub(p_var(ax),p_var(bx));
    auto dy=p_sub(p_var(ay),p_var(by));
    return f_atom(
        PolynomialRelation::GreaterZero,
        p_add({p_square(std::move(dx)),p_square(std::move(dy))}));
}

void collect_formula_stats(
    const SemialgebraicFormula& f,
    std::size_t& max_degree,
    std::size_t& atom_count) {
    if (f.kind==FormulaKind::Atom) {
        ++atom_count;
        max_degree=std::max(max_degree,polynomial_degree(f.atom.polynomial));
    }
    for (const auto& c:f.children)
        collect_formula_stats(c,max_degree,atom_count);
}

bool metadata_ok(const ExactSemialgebraicProgram& p) {
    return p.schema==kSchema &&
           p.schema_version==kSchemaVersion &&
           p.lowering_contract==kLoweringContract &&
           p.theorem_id==kTheoremId &&
           p.source_digest.algorithm=="SHA-256" &&
           p.source_digest.value==kSourceDigest;
}

std::string variable_scope_text(SemialgebraicVariableScope s) {
    return s==SemialgebraicVariableScope::Outer ? "outer" : "pair_local";
}

const char* algebraic_kind_text(ExactAlgebraicKind k) {
    switch (k) {
        case ExactAlgebraicKind::Rational: return "rational";
        case ExactAlgebraicKind::SinPi: return "sin_pi";
        case ExactAlgebraicKind::CosPi: return "cos_pi";
        case ExactAlgebraicKind::Add: return "add";
        case ExactAlgebraicKind::Multiply: return "mul";
        case ExactAlgebraicKind::Negate: return "neg";
        case ExactAlgebraicKind::Inverse: return "inv";
    }
    return "invalid";
}

const char* polynomial_kind_text(PolynomialExprKind k) {
    switch (k) {
        case PolynomialExprKind::Constant: return "const";
        case PolynomialExprKind::Variable: return "var";
        case PolynomialExprKind::Add: return "add";
        case PolynomialExprKind::Multiply: return "mul";
        case PolynomialExprKind::Negate: return "neg";
    }
    return "invalid";
}

const char* relation_text(PolynomialRelation r) {
    switch (r) {
        case PolynomialRelation::EqualZero: return "eq0";
        case PolynomialRelation::GreaterEqualZero: return "ge0";
        case PolynomialRelation::GreaterZero: return "gt0";
    }
    return "invalid";
}

const char* formula_kind_text(FormulaKind k) {
    switch (k) {
        case FormulaKind::True: return "true";
        case FormulaKind::False: return "false";
        case FormulaKind::Atom: return "atom";
        case FormulaKind::And: return "and";
        case FormulaKind::Or: return "or";
        case FormulaKind::Not: return "not";
        case FormulaKind::Exists: return "exists";
    }
    return "invalid";
}

void append_size_vector(std::ostringstream& os,
                        const std::vector<std::size_t>& xs) {
    os << "[";
    for (std::size_t i=0;i<xs.size();++i) {
        if (i) os << ",";
        os << xs[i];
    }
    os << "]";
}

void append_input(std::ostringstream& os,const NetworkClosureInput& in) {
    os << "{\"vertex_count\":" << in.vertex_count
       << ",\"trace_vertices\":";
    append_size_vector(os,in.trace_vertices);
    os << ",\"turns_pi\":[";
    for (std::size_t i=0;i<in.turns.size();++i) {
        if (i) os << ",";
        os << "\"" << rational_string(in.turns[i].pi) << "\"";
    }
    os << "]}";
}

} // namespace

std::size_t polynomial_degree(const PolynomialExpr& p) {
    switch (p.kind) {
        case PolynomialExprKind::Constant:
            return 0;
        case PolynomialExprKind::Variable:
            return 1;
        case PolynomialExprKind::Add: {
            std::size_t d=0;
            for (const auto& x:p.args) d=std::max(d,polynomial_degree(x));
            return d;
        }
        case PolynomialExprKind::Multiply: {
            std::size_t d=0;
            for (const auto& x:p.args) d+=polynomial_degree(x);
            return d;
        }
        case PolynomialExprKind::Negate:
            if (p.args.size()!=1)
                throw std::invalid_argument("negated polynomial must have one child");
            return polynomial_degree(p.args.front());
    }
    throw std::invalid_argument("unknown polynomial expression kind");
}

ExactSemialgebraicProgram lower_fixed_turn_semialgebraic_problem(
    const FixedTurnSemialgebraicProblem& ir) {

    const auto iv=verify_fixed_turn_semialgebraic_problem(ir);
    if (iv.status!=SemialgebraicProblemValidationStatus::Valid)
        throw std::invalid_argument("cannot lower invalid Phase 6B2A IR");

    ExactSemialgebraicProgram out;
    out.input=ir.input;
    out.source_ir_digest=fixed_turn_semialgebraic_problem_digest(ir);

    const std::size_t n=ir.input.vertex_count;
    const std::size_t m=ir.input.turns.size();

    out.variables.reserve(ir.outer_real_variable_count+
                          ir.pair_local_real_variable_count);
    out.outer_variables.reserve(ir.outer_real_variable_count);

    for (std::size_t v=0;v<n;++v) {
        const std::size_t ix=out.variables.size();
        out.variables.push_back({ix,"p"+std::to_string(v)+"_x",
                                 SemialgebraicVariableScope::Outer,0,0});
        out.outer_variables.push_back(ix);
        const std::size_t iy=out.variables.size();
        out.variables.push_back({iy,"p"+std::to_string(v)+"_y",
                                 SemialgebraicVariableScope::Outer,0,0});
        out.outer_variables.push_back(iy);
    }
    for (std::size_t e=0;e<m;++e) {
        const std::size_t i=out.variables.size();
        out.variables.push_back({i,"c"+std::to_string(e),
                                 SemialgebraicVariableScope::Outer,0,0});
        out.outer_variables.push_back(i);
    }

    std::vector<EdgePolynomials> derived;
    derived.reserve(m);
    for (const auto& e:ir.edges)
        derived.push_back(derive_edge(e,n));

    std::vector<SemialgebraicFormula> body;

    // Translation normalization p_root = 0.
    body.push_back(f_atom(
        PolynomialRelation::EqualZero,p_var(px(ir.root_vertex))));
    body.push_back(f_atom(
        PolynomialRelation::EqualZero,p_var(py(ir.root_vertex))));

    // Network equations p_head-p_tail=c*u, with exact algebraic u.
    for (const auto& e:ir.edges) {
        const auto ca=a_cos(e.chord_phase_pi);
        const auto sa=a_sin(e.chord_phase_pi);
        const auto c=p_var(cvar(n,e.edge));

        body.push_back(f_atom(
            PolynomialRelation::EqualZero,
            p_sub(
                p_sub(p_var(px(e.head)),p_var(px(e.tail))),
                p_scale(ca,c))));
        body.push_back(f_atom(
            PolynomialRelation::EqualZero,
            p_sub(
                p_sub(p_var(py(e.head)),p_var(py(e.tail))),
                p_scale(sa,c))));
    }

    // Positive metric and scale normalization.
    std::vector<PolynomialExpr> sum_c;
    sum_c.reserve(m+1);
    for (std::size_t e=0;e<m;++e) {
        body.push_back(f_atom(
            PolynomialRelation::GreaterZero,p_var(cvar(n,e))));
        sum_c.push_back(p_var(cvar(n,e)));
    }
    sum_c.push_back(p_int(-1));
    body.push_back(f_atom(
        PolynomialRelation::EqualZero,p_add(std::move(sum_c))));

    // Distinct quotient vertices.
    for (const auto& [v,w]:ir.distinct_vertex_pairs) {
        body.push_back(distinct_point(px(v),py(v),px(w),py(w)));
    }

    // Prescribed crossing transversality is fixed by the turn data. A failed
    // check makes the theorem sentence false before any metric search.
    for (const auto& t:ir.prescribed_transversality)
        if (!t.transverse) body.push_back(f_false());

    // Every pair contributes !(exists q I_ef).
    for (const auto& pair:ir.forbidden_pairs) {
        const std::size_t qx=out.variables.size();
        out.variables.push_back({
            qx,
            "q_"+std::to_string(pair.e)+"_"+std::to_string(pair.f)+"_x",
            SemialgebraicVariableScope::PairLocal,pair.e,pair.f});
        const std::size_t qy=out.variables.size();
        out.variables.push_back({
            qy,
            "q_"+std::to_string(pair.e)+"_"+std::to_string(pair.f)+"_y",
            SemialgebraicVariableScope::PairLocal,pair.e,pair.f});

        std::vector<SemialgebraicFormula> incidence;
        incidence.push_back(arc_membership(
            ir.edges[pair.e],derived[pair.e],qx,qy,n));
        incidence.push_back(arc_membership(
            ir.edges[pair.f],derived[pair.f],qx,qy,n));

        for (const auto v:pair.allowed_common_vertices)
            incidence.push_back(distinct_point(qx,qy,px(v),py(v)));

        body.push_back(f_not(f_exists(
            {qx,qy},
            f_and(std::move(incidence)))));
    }

    if (out.variables.size()!=
        ir.outer_real_variable_count+ir.pair_local_real_variable_count)
        throw std::runtime_error("Phase 6B2B variable accounting mismatch");

    out.formula=f_exists(
        out.outer_variables,
        f_and(std::move(body)));

    collect_formula_stats(
        out.formula,out.maximum_polynomial_degree,out.polynomial_atom_count);

    if (out.maximum_polynomial_degree>2)
        throw std::runtime_error(
            "theorem lowering unexpectedly exceeded quadratic degree");

    return out;
}

ExactSemialgebraicProgram compile_exact_semialgebraic_program(
    const NetworkClosureInput& input) {
    return lower_fixed_turn_semialgebraic_problem(
        compile_fixed_turn_semialgebraic_problem(input));
}

std::string serialize_exact_algebraic_expr(
    const ExactAlgebraicExpr& x) {
    std::ostringstream os;
    os << "{\"kind\":\"" << algebraic_kind_text(x.kind) << "\"";
    if (x.kind==ExactAlgebraicKind::Rational)
        os << ",\"value\":\"" << rational_text(x.rational) << "\"";
    if (x.kind==ExactAlgebraicKind::SinPi ||
        x.kind==ExactAlgebraicKind::CosPi)
        os << ",\"phase_pi\":\"" << rational_string(x.phase_pi) << "\"";
    if (!x.args.empty()) {
        os << ",\"args\":[";
        for (std::size_t i=0;i<x.args.size();++i) {
            if (i) os << ",";
            os << serialize_exact_algebraic_expr(x.args[i]);
        }
        os << "]";
    }
    os << "}";
    return os.str();
}

std::string serialize_polynomial_expr(
    const PolynomialExpr& p) {
    std::ostringstream os;
    os << "{\"kind\":\"" << polynomial_kind_text(p.kind) << "\"";
    if (p.kind==PolynomialExprKind::Constant)
        os << ",\"value\":" << serialize_exact_algebraic_expr(p.constant);
    if (p.kind==PolynomialExprKind::Variable)
        os << ",\"variable\":" << p.variable;
    if (!p.args.empty()) {
        os << ",\"args\":[";
        for (std::size_t i=0;i<p.args.size();++i) {
            if (i) os << ",";
            os << serialize_polynomial_expr(p.args[i]);
        }
        os << "]";
    }
    os << "}";
    return os.str();
}

std::string serialize_semialgebraic_formula(
    const SemialgebraicFormula& f) {
    std::ostringstream os;
    os << "{\"kind\":\"" << formula_kind_text(f.kind) << "\"";
    if (f.kind==FormulaKind::Atom) {
        os << ",\"relation\":\"" << relation_text(f.atom.relation)
           << "\",\"polynomial\":"
           << serialize_polynomial_expr(f.atom.polynomial);
    }
    if (f.kind==FormulaKind::Exists) {
        os << ",\"variables\":";
        append_size_vector(os,f.quantified_variables);
    }
    if (!f.children.empty()) {
        os << ",\"children\":[";
        for (std::size_t i=0;i<f.children.size();++i) {
            if (i) os << ",";
            os << serialize_semialgebraic_formula(f.children[i]);
        }
        os << "]";
    }
    os << "}";
    return os.str();
}

std::string serialize_exact_semialgebraic_program(
    const ExactSemialgebraicProgram& p) {
    std::ostringstream os;
    os << "{\"schema\":\"" << p.schema
       << "\",\"schema_version\":\"" << p.schema_version
       << "\",\"lowering_contract\":\"" << p.lowering_contract
       << "\",\"theorem_id\":\"" << p.theorem_id
       << "\",\"source_digest\":{\"algorithm\":\"" << p.source_digest.algorithm
       << "\",\"value\":\"" << p.source_digest.value
       << "\"},\"source_ir_digest\":\"" << p.source_ir_digest
       << "\",\"input\":";
    append_input(os,p.input);

    os << ",\"variables\":[";
    for (std::size_t i=0;i<p.variables.size();++i) {
        if (i) os << ",";
        const auto& v=p.variables[i];
        os << "{\"index\":" << v.index
           << ",\"name\":\"" << v.name
           << "\",\"scope\":\"" << variable_scope_text(v.scope) << "\"";
        if (v.scope==SemialgebraicVariableScope::PairLocal)
            os << ",\"pair\":[" << v.pair_e << "," << v.pair_f << "]";
        os << "}";
    }
    os << "],\"outer_variables\":";
    append_size_vector(os,p.outer_variables);
    os << ",\"maximum_polynomial_degree\":" << p.maximum_polynomial_degree
       << ",\"polynomial_atom_count\":" << p.polynomial_atom_count
       << ",\"formula\":" << serialize_semialgebraic_formula(p.formula)
       << "}";
    return os.str();
}

std::string exact_semialgebraic_program_digest(
    const ExactSemialgebraicProgram& program) {
    return sha256_hex(serialize_exact_semialgebraic_program(program));
}

ExactSemialgebraicProgramValidation verify_exact_semialgebraic_program(
    const ExactSemialgebraicProgram& program) {
    try {
        if (!metadata_ok(program))
            return {
                ExactSemialgebraicProgramValidationStatus::InvalidMetadata,
                "schema, lowering contract, theorem id, or source digest mismatch"};

        const auto ir=compile_fixed_turn_semialgebraic_problem(program.input);
        const auto ir_validation=verify_fixed_turn_semialgebraic_problem(ir);
        if (ir_validation.status!=SemialgebraicProblemValidationStatus::Valid)
            return {
                ExactSemialgebraicProgramValidationStatus::InvalidSourceIrBinding,
                "recompiled Phase 6B2A IR is invalid"};

        if (program.source_ir_digest!=fixed_turn_semialgebraic_problem_digest(ir))
            return {
                ExactSemialgebraicProgramValidationStatus::InvalidSourceIrBinding,
                "Phase 6B2A IR digest mismatch"};

        std::size_t max_degree=0, atom_count=0;
        collect_formula_stats(program.formula,max_degree,atom_count);
        if (max_degree>2 || program.maximum_polynomial_degree!=max_degree)
            return {
                ExactSemialgebraicProgramValidationStatus::DegreeContractViolation,
                "polynomial degree exceeds or disagrees with the quadratic theorem contract"};
        if (program.polynomial_atom_count!=atom_count)
            return {
                ExactSemialgebraicProgramValidationStatus::InvalidStructure,
                "polynomial atom count mismatch"};

        const auto expected=lower_fixed_turn_semialgebraic_problem(ir);
        if (serialize_exact_semialgebraic_program(program)!=
            serialize_exact_semialgebraic_program(expected))
            return {
                ExactSemialgebraicProgramValidationStatus::InvalidStructure,
                "program differs from deterministic exact theorem lowering"};

        return {ExactSemialgebraicProgramValidationStatus::Valid,""};
    } catch (const std::exception& e) {
        return {
            ExactSemialgebraicProgramValidationStatus::InvalidStructure,
            e.what()};
    } catch (...) {
        return {
            ExactSemialgebraicProgramValidationStatus::InvalidStructure,
            "unknown exact semialgebraic program verification failure"};
    }
}

} // namespace pcg
