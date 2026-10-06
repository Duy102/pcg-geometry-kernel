#include <pcg/qe_backend.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pcg {
namespace {

constexpr const char* kSchema="pcg-qepcad-rcf-program";
constexpr const char* kSchemaVersion="1.0";
constexpr const char* kBackendContract="PCG-FT-QE-QEPCAD-001";
constexpr const char* kTheoremId="PCG-GENERAL-FIXED-TURN-INTERSECTION-THM-001";
constexpr const char* kSourceDigest=
    "625f91fac3d201f29b7a8a253b0e8fb3de5a8358f5eaeb2527ee8a85899ce208";
constexpr unsigned kMaximumTrigDenominator=256;

using RatPoly=std::vector<Rational>;

std::string rational_text(const Rational& r) {
    std::ostringstream os;
    os << r.numerator() << "/" << r.denominator();
    return os.str();
}

std::string json_escape(const std::string& s) {
    std::ostringstream os;
    for (char ch:s) {
        switch (ch) {
            case '\\': os << "\\\\"; break;
            case '"': os << "\\\""; break;
            case '\n': os << "\\n"; break;
            case '\r': os << "\\r"; break;
            case '\t': os << "\\t"; break;
            default: os << ch; break;
        }
    }
    return os.str();
}

void trim(RatPoly& p) {
    while (p.size()>1 && p.back()==Rational(BigInt{0}))
        p.pop_back();
    if (p.empty()) p.push_back(Rational(BigInt{0}));
}

bool zero_poly(const RatPoly& p) {
    for (const auto& x:p)
        if (x!=Rational(BigInt{0})) return false;
    return true;
}

RatPoly derivative(const RatPoly& p) {
    if (p.size()<=1) return {Rational(BigInt{0})};
    RatPoly d(p.size()-1,Rational(BigInt{0}));
    for (std::size_t i=1;i<p.size();++i)
        d[i-1]=p[i]*Rational(BigInt{i});
    trim(d);
    return d;
}

std::pair<RatPoly,RatPoly> divrem(RatPoly a,RatPoly b) {
    trim(a); trim(b);
    if (zero_poly(b))
        throw std::invalid_argument("zero polynomial divisor");
    RatPoly q(
        a.size()>=b.size()?a.size()-b.size()+1:1,
        Rational(BigInt{0}));
    while (!zero_poly(a) && a.size()>=b.size()) {
        const std::size_t shift=a.size()-b.size();
        const Rational c=a.back()/b.back();
        q[shift]+=c;
        for (std::size_t j=0;j<b.size();++j)
            a[shift+j]-=c*b[j];
        trim(a);
    }
    trim(q); trim(a);
    return {q,a};
}

RatPoly monic(RatPoly p) {
    trim(p);
    if (zero_poly(p)) return p;
    const Rational lead=p.back();
    for (auto& x:p) x/=lead;
    return p;
}

RatPoly gcd_poly(RatPoly a,RatPoly b) {
    trim(a); trim(b);
    while (!zero_poly(b)) {
        auto qr=divrem(a,b);
        a=std::move(b);
        b=std::move(qr.second);
    }
    return monic(std::move(a));
}

RatPoly exact_quotient(RatPoly a,const RatPoly& b) {
    auto qr=divrem(std::move(a),b);
    if (!zero_poly(qr.second))
        throw std::runtime_error("non-exact polynomial quotient");
    return qr.first;
}

RatPoly square_free_part(RatPoly p) {
    trim(p);
    const RatPoly d=derivative(p);
    if (zero_poly(d)) return monic(std::move(p));
    const RatPoly g=gcd_poly(p,d);
    return monic(exact_quotient(std::move(p),g));
}

Rational eval_poly(const RatPoly& p,const Rational& x) {
    Rational y(BigInt{0});
    for (std::size_t i=p.size();i-- >0;)
        y=y*x+p[i];
    return y;
}

int rational_sign(const Rational& x) {
    if (x<Rational(BigInt{0})) return -1;
    if (x>Rational(BigInt{0})) return +1;
    return 0;
}

std::vector<RatPoly> sturm_sequence(RatPoly p) {
    p=square_free_part(std::move(p));
    if (zero_poly(p))
        throw std::invalid_argument("zero polynomial has no finite root isolation");
    std::vector<RatPoly> seq;
    seq.push_back(p);
    RatPoly d=derivative(p);
    if (zero_poly(d)) return seq;
    seq.push_back(std::move(d));
    while (true) {
        auto qr=divrem(seq[seq.size()-2],seq.back());
        RatPoly r=std::move(qr.second);
        if (zero_poly(r)) break;
        for (auto& x:r) x=-x;
        trim(r);
        seq.push_back(std::move(r));
    }
    return seq;
}

int sturm_variations(
    const std::vector<RatPoly>& seq,
    const Rational& x) {
    int previous=0;
    int variations=0;
    for (const auto& p:seq) {
        const int s=rational_sign(eval_poly(p,x));
        if (s==0) continue;
        if (previous!=0 && s!=previous) ++variations;
        previous=s;
    }
    return variations;
}

std::size_t root_count_open(
    const IntegerPolynomial& p,
    const Rational& lo,
    const Rational& hi) {
    if (!(lo<hi))
        throw std::invalid_argument("root isolation interval must be ordered");
    RatPoly q;
    q.reserve(p.coefficients.size());
    for (const auto& x:p.coefficients)
        q.push_back(Rational(x));
    trim(q);
    RatPoly sf=square_free_part(q);
    if (eval_poly(sf,lo)==Rational(BigInt{0}) ||
        eval_poly(sf,hi)==Rational(BigInt{0}))
        throw std::invalid_argument("Sturm interval endpoint is a root");
    const auto seq=sturm_sequence(std::move(q));
    const int n=sturm_variations(seq,lo)-sturm_variations(seq,hi);
    if (n<0)
        throw std::runtime_error("negative Sturm root count");
    return static_cast<std::size_t>(n);
}

std::vector<BigInt> poly_sub(
    std::vector<BigInt> a,
    const std::vector<BigInt>& b) {
    if (a.size()<b.size()) a.resize(b.size(),BigInt{0});
    for (std::size_t i=0;i<b.size();++i) a[i]-=b[i];
    while (a.size()>1 && a.back()==0) a.pop_back();
    return a;
}

std::vector<BigInt> poly_scale(
    std::vector<BigInt> a,
    const BigInt& s) {
    for (auto& x:a) x*=s;
    while (a.size()>1 && a.back()==0) a.pop_back();
    return a;
}

std::vector<BigInt> poly_shift_x(std::vector<BigInt> a) {
    a.insert(a.begin(),BigInt{0});
    return a;
}

IntegerPolynomial chebyshev_relation(PiRational phase) {
    const BigInt den_big=phase.value.denominator();
    if (den_big>BigInt{kMaximumTrigDenominator})
        throw std::invalid_argument(
            "QEPCAD adapter trig denominator exceeds exact isolation limit");
    const unsigned d=den_big.convert_to<unsigned>();
    if (d==0) throw std::invalid_argument("zero trig denominator");

    std::vector<BigInt> t0{BigInt{1}};
    if (d==0) return IntegerPolynomial{t0};
    std::vector<BigInt> t1{BigInt{0},BigInt{1}};
    std::vector<BigInt> td;
    if (d==1) td=t1;
    else {
        for (unsigned k=2;k<=d;++k) {
            auto next=poly_sub(
                poly_scale(poly_shift_x(t1),BigInt{2}),t0);
            t0=std::move(t1);
            t1=std::move(next);
        }
        td=t1;
    }

    BigInt rem=phase.value.numerator()%2;
    if (rem<0) rem=-rem;
    const BigInt target=(rem==0)?BigInt{1}:BigInt{-1};
    td[0]-=target;
    while (td.size()>1 && td.back()==0) td.pop_back();
    return IntegerPolynomial{std::move(td)};
}

unsigned target_root_rank(PiRational phase) {
    const BigInt d_big=phase.value.denominator();
    if (d_big>BigInt{kMaximumTrigDenominator})
        throw std::invalid_argument("trig denominator exceeds root-rank limit");
    const unsigned d=d_big.convert_to<unsigned>();
    const BigInt period=BigInt{2}*d_big;
    BigInt n=phase.value.numerator()%period;
    if (n<0) n+=period;
    BigInt reflected=period-n;
    BigInt k_big=n<reflected?n:reflected;
    const unsigned k=k_big.convert_to<unsigned>();

    unsigned jmax=d;
    if (((jmax-k)&1U)!=0U) {
        if (jmax==0)
            throw std::runtime_error("invalid cosine root parity");
        --jmax;
    }
    if (jmax<k)
        throw std::runtime_error("cosine root rank underflow");
    return (jmax-k)/2U+1U;
}

std::pair<Rational,Rational> isolate_root_by_rank(
    const IntegerPolynomial& p,
    unsigned rank) {
    Rational lo(BigInt{-2});
    Rational hi(BigInt{2});
    std::size_t total=root_count_open(p,lo,hi);
    if (rank==0 || rank>total)
        throw std::runtime_error("target algebraic root rank is out of range");

    std::size_t local_rank=rank;
    for (unsigned iteration=0;iteration<2048;++iteration) {
        const std::size_t current=root_count_open(p,lo,hi);
        if (current==1) return {lo,hi};
        if (local_rank==0 || local_rank>current)
            throw std::runtime_error("root rank escaped isolation interval");

        Rational split=(lo+hi)/Rational(BigInt{2});

        RatPoly q;
        q.reserve(p.coefficients.size());
        for (const auto& x:p.coefficients) q.push_back(Rational(x));
        RatPoly sf=square_free_part(std::move(q));

        if (eval_poly(sf,split)==Rational(BigInt{0})) {
            bool found=false;
            for (unsigned m=2;m<=64;++m) {
                split=(lo*Rational(BigInt{m})+hi)/
                      Rational(BigInt{m+1});
                if (eval_poly(sf,split)!=Rational(BigInt{0})) {
                    found=true;
                    break;
                }
            }
            if (!found)
                throw std::runtime_error(
                    "failed to choose a non-root rational split point");
        }

        const std::size_t left=root_count_open(p,lo,split);
        if (local_rank<=left) {
            hi=split;
        } else {
            local_rank-=left;
            lo=split;
        }
    }
    throw std::runtime_error("exact algebraic root isolation iteration limit");
}

TrigIsolationDefinition build_trig_definition(
    PiRational phase,
    std::size_t cosine_variable,
    std::size_t sine_variable) {
    IntegerPolynomial p=chebyshev_relation(phase);
    const unsigned rank=target_root_rank(phase);
    auto interval=isolate_root_by_rank(p,rank);
    if (root_count_open(p,interval.first,interval.second)!=1)
        throw std::runtime_error("cosine isolation did not certify one root");
    return {
        phase,
        cosine_variable,
        sine_variable,
        std::move(p),
        std::move(interval.first),
        std::move(interval.second),
        sin_pi_sign(phase)
    };
}

RcfPolynomialExpr r_rat(Rational q) {
    RcfPolynomialExpr x;
    x.kind=RcfPolynomialKind::Rational;
    x.rational=std::move(q);
    return x;
}

RcfPolynomialExpr r_int(long long n) {
    return r_rat(Rational(BigInt{n}));
}

RcfPolynomialExpr r_var(std::size_t v) {
    RcfPolynomialExpr x;
    x.kind=RcfPolynomialKind::Variable;
    x.variable=v;
    return x;
}

RcfPolynomialExpr r_add(std::vector<RcfPolynomialExpr> xs) {
    if (xs.empty()) return r_int(0);
    if (xs.size()==1) return std::move(xs.front());
    RcfPolynomialExpr x;
    x.kind=RcfPolynomialKind::Add;
    x.args=std::move(xs);
    return x;
}

RcfPolynomialExpr r_mul(std::vector<RcfPolynomialExpr> xs) {
    if (xs.empty()) return r_int(1);
    if (xs.size()==1) return std::move(xs.front());
    RcfPolynomialExpr x;
    x.kind=RcfPolynomialKind::Multiply;
    x.args=std::move(xs);
    return x;
}

RcfPolynomialExpr r_neg(RcfPolynomialExpr a) {
    RcfPolynomialExpr x;
    x.kind=RcfPolynomialKind::Negate;
    x.args.push_back(std::move(a));
    return x;
}

RcfPolynomialExpr r_sub(RcfPolynomialExpr a,RcfPolynomialExpr b) {
    return r_add({std::move(a),r_neg(std::move(b))});
}

RcfPolynomialExpr r_square(RcfPolynomialExpr a) {
    RcfPolynomialExpr b=a;
    return r_mul({std::move(a),std::move(b)});
}

RcfFormula rf_true() {
    RcfFormula f;
    f.kind=RcfFormulaKind::True;
    return f;
}

RcfFormula rf_false() {
    RcfFormula f;
    f.kind=RcfFormulaKind::False;
    return f;
}

RcfFormula rf_atom(RcfRelation rel,RcfPolynomialExpr p) {
    RcfFormula f;
    f.kind=RcfFormulaKind::Atom;
    f.atom.relation=rel;
    f.atom.polynomial=std::move(p);
    return f;
}

RcfFormula rf_and(std::vector<RcfFormula> xs) {
    if (xs.empty()) return rf_true();
    if (xs.size()==1) return std::move(xs.front());
    RcfFormula f;
    f.kind=RcfFormulaKind::And;
    f.children=std::move(xs);
    return f;
}

RcfFormula rf_or(std::vector<RcfFormula> xs) {
    if (xs.empty()) return rf_false();
    if (xs.size()==1) return std::move(xs.front());
    RcfFormula f;
    f.kind=RcfFormulaKind::Or;
    f.children=std::move(xs);
    return f;
}

RcfFormula rf_not(RcfFormula a) {
    RcfFormula f;
    f.kind=RcfFormulaKind::Not;
    f.children.push_back(std::move(a));
    return f;
}

RcfFormula rf_exists(
    std::vector<std::size_t> vars,
    RcfFormula body) {
    if (vars.empty()) return body;
    RcfFormula f;
    f.kind=RcfFormulaKind::Exists;
    f.quantified_variables=std::move(vars);
    f.children.push_back(std::move(body));
    return f;
}

RcfPolynomialExpr integer_polynomial_expr(
    const IntegerPolynomial& p,
    std::size_t variable) {
    RcfPolynomialExpr acc=r_int(0);
    for (std::size_t i=p.coefficients.size();i-- >0;) {
        acc=r_add({
            r_mul({std::move(acc),r_var(variable)}),
            r_rat(Rational(p.coefficients[i]))
        });
    }
    return acc;
}

const char* variable_role_text(QeVariableRole r) {
    switch (r) {
        case QeVariableRole::AlgebraicCosine: return "algebraic_cosine";
        case QeVariableRole::AlgebraicSine: return "algebraic_sine";
        case QeVariableRole::AlgebraicOperation: return "algebraic_operation";
        case QeVariableRole::SourceOuter: return "source_outer";
        case QeVariableRole::SourcePairLocal: return "source_pair_local";
    }
    return "invalid";
}

const char* poly_kind_text(RcfPolynomialKind k) {
    switch (k) {
        case RcfPolynomialKind::Rational: return "rational";
        case RcfPolynomialKind::Variable: return "variable";
        case RcfPolynomialKind::Add: return "add";
        case RcfPolynomialKind::Multiply: return "mul";
        case RcfPolynomialKind::Negate: return "neg";
    }
    return "invalid";
}

const char* relation_text(RcfRelation r) {
    switch (r) {
        case RcfRelation::EqualZero: return "eq0";
        case RcfRelation::GreaterEqualZero: return "ge0";
        case RcfRelation::GreaterZero: return "gt0";
    }
    return "invalid";
}

const char* formula_kind_text(RcfFormulaKind k) {
    switch (k) {
        case RcfFormulaKind::True: return "true";
        case RcfFormulaKind::False: return "false";
        case RcfFormulaKind::Atom: return "atom";
        case RcfFormulaKind::And: return "and";
        case RcfFormulaKind::Or: return "or";
        case RcfFormulaKind::Not: return "not";
        case RcfFormulaKind::Exists: return "exists";
    }
    return "invalid";
}

void append_size_vector(
    std::ostringstream& os,
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

class RcfLowerer {
public:
    explicit RcfLowerer(const ExactSemialgebraicProgram& source)
        : source_(source) {}

    QepcadRcfProgram run() {
        QepcadRcfProgram out;
        out.input=source_.input;
        out.source_program_digest=exact_semialgebraic_program_digest(source_);
        out_= &out;

        scan_formula(source_.formula);

        out.coefficient_variables.reserve(out.variables.size());
        for (std::size_t i=0;i<out.variables.size();++i)
            out.coefficient_variables.push_back(i);

        const std::size_t source_offset=out.variables.size();
        source_map_.resize(source_.variables.size());
        for (const auto& v:source_.variables) {
            const std::size_t idx=out.variables.size();
            source_map_[v.index]=idx;
            QeVariableRole role=
                v.scope==SemialgebraicVariableScope::Outer
                ?QeVariableRole::SourceOuter
                :QeVariableRole::SourcePairLocal;
            out.variables.push_back({
                idx,
                v.name,
                role,
                PiRational{0,1},
                v.index,
                v.pair_e,
                v.pair_f
            });
        }
        if (out.variables.size()!=source_offset+source_.variables.size())
            throw std::runtime_error("QE source variable accounting mismatch");

        std::vector<RcfFormula> definitions;
        for (const auto& t:out.trig_definitions) {
            definitions.push_back(rf_atom(
                RcfRelation::EqualZero,
                integer_polynomial_expr(
                    t.cosine_polynomial,t.cosine_variable)));
            definitions.push_back(rf_atom(
                RcfRelation::EqualZero,
                r_add({
                    r_square(r_var(t.cosine_variable)),
                    r_square(r_var(t.sine_variable)),
                    r_int(-1)
                })));
            definitions.push_back(rf_atom(
                RcfRelation::GreaterZero,
                r_sub(
                    r_var(t.cosine_variable),
                    r_rat(t.cosine_lower))));
            definitions.push_back(rf_atom(
                RcfRelation::GreaterZero,
                r_sub(
                    r_rat(t.cosine_upper),
                    r_var(t.cosine_variable))));
            if (t.sine_sign>0)
                definitions.push_back(rf_atom(
                    RcfRelation::GreaterZero,
                    r_var(t.sine_variable)));
            else if (t.sine_sign<0)
                definitions.push_back(rf_atom(
                    RcfRelation::GreaterZero,
                    r_neg(r_var(t.sine_variable))));
            else
                definitions.push_back(rf_atom(
                    RcfRelation::EqualZero,
                    r_var(t.sine_variable)));
        }

        for (const auto& node:operation_nodes_)
            definitions.push_back(operation_definition(node));

        definitions.push_back(translate_formula(source_.formula));
        out.formula=rf_exists(
            out.coefficient_variables,
            rf_and(std::move(definitions)));

        collect_stats(
            out.formula,out.maximum_total_degree,out.rational_atom_count);
        out_=nullptr;
        return out;
    }

private:
    struct OperationNode {
        std::size_t variable{};
        ExactAlgebraicExpr expression;
        std::string serialization;
    };

    const ExactSemialgebraicProgram& source_;
    QepcadRcfProgram* out_{};
    std::map<std::string,std::pair<std::size_t,std::size_t>> trig_;
    std::map<std::string,std::size_t> operations_;
    std::vector<OperationNode> operation_nodes_;
    std::vector<std::size_t> source_map_;

    void register_trig(PiRational phase) {
        const std::string key=rational_string(phase);
        if (trig_.find(key)!=trig_.end()) return;

        const std::size_t c=out_->variables.size();
        out_->variables.push_back({
            c,"cos_"+std::to_string(c),
            QeVariableRole::AlgebraicCosine,
            phase,0,0,0
        });
        const std::size_t s=out_->variables.size();
        out_->variables.push_back({
            s,"sin_"+std::to_string(s),
            QeVariableRole::AlgebraicSine,
            phase,0,0,0
        });
        trig_[key]={c,s};
        out_->trig_definitions.push_back(
            build_trig_definition(phase,c,s));
    }

    void scan_alg(const ExactAlgebraicExpr& a) {
        switch (a.kind) {
            case ExactAlgebraicKind::Rational:
                return;
            case ExactAlgebraicKind::SinPi:
            case ExactAlgebraicKind::CosPi:
                register_trig(a.phase_pi);
                return;
            case ExactAlgebraicKind::Add:
            case ExactAlgebraicKind::Multiply:
            case ExactAlgebraicKind::Negate:
            case ExactAlgebraicKind::Inverse:
                for (const auto& x:a.args) scan_alg(x);
                break;
        }

        const std::string key=serialize_exact_algebraic_expr(a);
        if (operations_.find(key)!=operations_.end()) return;
        const std::size_t idx=out_->variables.size();
        out_->variables.push_back({
            idx,
            "alg_"+std::to_string(idx),
            QeVariableRole::AlgebraicOperation,
            PiRational{0,1},0,0,0
        });
        operations_[key]=idx;
        operation_nodes_.push_back({idx,a,key});
        out_->operation_definitions.push_back({idx,key});
    }

    void scan_poly(const PolynomialExpr& p) {
        if (p.kind==PolynomialExprKind::Constant)
            scan_alg(p.constant);
        for (const auto& x:p.args) scan_poly(x);
    }

    void scan_formula(const SemialgebraicFormula& f) {
        if (f.kind==FormulaKind::Atom)
            scan_poly(f.atom.polynomial);
        for (const auto& c:f.children) scan_formula(c);
    }

    RcfPolynomialExpr alg_term(const ExactAlgebraicExpr& a) const {
        if (a.kind==ExactAlgebraicKind::Rational)
            return r_rat(a.rational);
        if (a.kind==ExactAlgebraicKind::SinPi ||
            a.kind==ExactAlgebraicKind::CosPi) {
            const auto it=trig_.find(rational_string(a.phase_pi));
            if (it==trig_.end())
                throw std::runtime_error("missing registered trig atom");
            return r_var(
                a.kind==ExactAlgebraicKind::CosPi
                ?it->second.first
                :it->second.second);
        }
        const auto it=operations_.find(
            serialize_exact_algebraic_expr(a));
        if (it==operations_.end())
            throw std::runtime_error("missing registered algebraic operation");
        return r_var(it->second);
    }

    RcfFormula operation_definition(const OperationNode& n) const {
        const auto& a=n.expression;
        const auto z=r_var(n.variable);
        switch (a.kind) {
            case ExactAlgebraicKind::Add: {
                std::vector<RcfPolynomialExpr> xs;
                for (const auto& x:a.args) xs.push_back(alg_term(x));
                return rf_atom(
                    RcfRelation::EqualZero,
                    r_sub(z,r_add(std::move(xs))));
            }
            case ExactAlgebraicKind::Multiply: {
                std::vector<RcfPolynomialExpr> xs;
                for (const auto& x:a.args) xs.push_back(alg_term(x));
                return rf_atom(
                    RcfRelation::EqualZero,
                    r_sub(z,r_mul(std::move(xs))));
            }
            case ExactAlgebraicKind::Negate:
                if (a.args.size()!=1)
                    throw std::invalid_argument(
                        "algebraic negate must have one argument");
                return rf_atom(
                    RcfRelation::EqualZero,
                    r_sub(z,r_neg(alg_term(a.args.front()))));
            case ExactAlgebraicKind::Inverse:
                if (a.args.size()!=1)
                    throw std::invalid_argument(
                        "algebraic inverse must have one argument");
                return rf_atom(
                    RcfRelation::EqualZero,
                    r_sub(
                        r_mul({z,alg_term(a.args.front())}),
                        r_int(1)));
            case ExactAlgebraicKind::Rational:
            case ExactAlgebraicKind::SinPi:
            case ExactAlgebraicKind::CosPi:
                break;
        }
        throw std::invalid_argument(
            "unexpected atomic algebraic operation definition");
    }

    RcfPolynomialExpr translate_poly(const PolynomialExpr& p) const {
        switch (p.kind) {
            case PolynomialExprKind::Constant:
                return alg_term(p.constant);
            case PolynomialExprKind::Variable:
                if (p.variable>=source_map_.size())
                    throw std::invalid_argument(
                        "source polynomial variable index out of range");
                return r_var(source_map_[p.variable]);
            case PolynomialExprKind::Add: {
                std::vector<RcfPolynomialExpr> xs;
                for (const auto& x:p.args)
                    xs.push_back(translate_poly(x));
                return r_add(std::move(xs));
            }
            case PolynomialExprKind::Multiply: {
                std::vector<RcfPolynomialExpr> xs;
                for (const auto& x:p.args)
                    xs.push_back(translate_poly(x));
                return r_mul(std::move(xs));
            }
            case PolynomialExprKind::Negate:
                if (p.args.size()!=1)
                    throw std::invalid_argument(
                        "source polynomial negate must have one child");
                return r_neg(translate_poly(p.args.front()));
        }
        throw std::invalid_argument("unknown source polynomial kind");
    }

    RcfFormula translate_formula(const SemialgebraicFormula& f) const {
        switch (f.kind) {
            case FormulaKind::True:
                return rf_true();
            case FormulaKind::False:
                return rf_false();
            case FormulaKind::Atom: {
                RcfRelation rel=RcfRelation::EqualZero;
                if (f.atom.relation==PolynomialRelation::GreaterEqualZero)
                    rel=RcfRelation::GreaterEqualZero;
                else if (f.atom.relation==PolynomialRelation::GreaterZero)
                    rel=RcfRelation::GreaterZero;
                return rf_atom(
                    rel,translate_poly(f.atom.polynomial));
            }
            case FormulaKind::And:
            case FormulaKind::Or: {
                std::vector<RcfFormula> xs;
                for (const auto& c:f.children)
                    xs.push_back(translate_formula(c));
                return f.kind==FormulaKind::And
                    ?rf_and(std::move(xs))
                    :rf_or(std::move(xs));
            }
            case FormulaKind::Not:
                if (f.children.size()!=1)
                    throw std::invalid_argument(
                        "source formula not must have one child");
                return rf_not(
                    translate_formula(f.children.front()));
            case FormulaKind::Exists: {
                if (f.children.size()!=1)
                    throw std::invalid_argument(
                        "source formula exists must have one body");
                std::vector<std::size_t> vars;
                vars.reserve(f.quantified_variables.size());
                for (const auto v:f.quantified_variables) {
                    if (v>=source_map_.size())
                        throw std::invalid_argument(
                            "source quantified variable out of range");
                    vars.push_back(source_map_[v]);
                }
                return rf_exists(
                    std::move(vars),
                    translate_formula(f.children.front()));
            }
        }
        throw std::invalid_argument("unknown source formula kind");
    }

    static void collect_stats(
        const RcfFormula& f,
        std::size_t& max_degree,
        std::size_t& atom_count) {
        if (f.kind==RcfFormulaKind::Atom) {
            ++atom_count;
            max_degree=std::max(
                max_degree,
                rcf_polynomial_degree(f.atom.polynomial));
        }
        for (const auto& c:f.children)
            collect_stats(c,max_degree,atom_count);
    }
};

bool metadata_ok(const QepcadRcfProgram& p) {
    return p.schema==kSchema &&
           p.schema_version==kSchemaVersion &&
           p.backend_contract==kBackendContract &&
           p.theorem_id==kTheoremId &&
           p.source_digest.algorithm=="SHA-256" &&
           p.source_digest.value==kSourceDigest;
}

struct PrenexQuantifier {
    bool existential{};
    std::size_t variable{};
};

struct PrenexResult {
    std::vector<PrenexQuantifier> quantifiers;
    RcfFormula matrix;
};

PrenexResult prenex(const RcfFormula& f,bool negate=false) {
    switch (f.kind) {
        case RcfFormulaKind::True:
            return {{},negate?rf_false():rf_true()};
        case RcfFormulaKind::False:
            return {{},negate?rf_true():rf_false()};
        case RcfFormulaKind::Atom:
            return {{},negate?rf_not(f):f};
        case RcfFormulaKind::Not:
            if (f.children.size()!=1)
                throw std::invalid_argument("RCF not must have one child");
            return prenex(f.children.front(),!negate);
        case RcfFormulaKind::Exists: {
            if (f.children.size()!=1)
                throw std::invalid_argument("RCF exists must have one child");
            auto sub=prenex(f.children.front(),negate);
            std::vector<PrenexQuantifier> q;
            q.reserve(f.quantified_variables.size()+sub.quantifiers.size());
            for (const auto v:f.quantified_variables)
                q.push_back({!negate,v});
            q.insert(q.end(),sub.quantifiers.begin(),sub.quantifiers.end());
            sub.quantifiers=std::move(q);
            return sub;
        }
        case RcfFormulaKind::And:
        case RcfFormulaKind::Or: {
            const bool use_and=
                (f.kind==RcfFormulaKind::And) != negate;
            std::vector<PrenexQuantifier> qs;
            std::vector<RcfFormula> bodies;
            for (const auto& c:f.children) {
                auto sub=prenex(c,negate);
                qs.insert(
                    qs.end(),sub.quantifiers.begin(),sub.quantifiers.end());
                bodies.push_back(std::move(sub.matrix));
            }
            return {
                std::move(qs),
                use_and?rf_and(std::move(bodies)):rf_or(std::move(bodies))
            };
        }
    }
    throw std::invalid_argument("unknown RCF formula kind");
}

std::string qepcad_rational(const Rational& r) {
    std::ostringstream os;
    if (r.denominator()==BigInt{1})
        os << r.numerator();
    else
        os << "(" << r.numerator() << "/" << r.denominator() << ")";
    return os.str();
}

std::string qepcad_poly(const RcfPolynomialExpr& p) {
    switch (p.kind) {
        case RcfPolynomialKind::Rational:
            return qepcad_rational(p.rational);
        case RcfPolynomialKind::Variable:
            return "v"+std::to_string(p.variable);
        case RcfPolynomialKind::Negate:
            if (p.args.size()!=1)
                throw std::invalid_argument(
                    "RCF polynomial negation must have one child");
            return "(-1 "+qepcad_poly(p.args.front())+")";
        case RcfPolynomialKind::Add: {
            if (p.args.empty()) return "0";
            std::ostringstream os;
            os << "(";
            for (std::size_t i=0;i<p.args.size();++i) {
                if (i) os << " + ";
                os << qepcad_poly(p.args[i]);
            }
            os << ")";
            return os.str();
        }
        case RcfPolynomialKind::Multiply: {
            if (p.args.empty()) return "1";
            std::ostringstream os;
            os << "(";
            for (std::size_t i=0;i<p.args.size();++i) {
                if (i) os << " ";
                os << qepcad_poly(p.args[i]);
            }
            os << ")";
            return os.str();
        }
    }
    throw std::invalid_argument("unknown RCF polynomial kind");
}

std::string qepcad_matrix(const RcfFormula& f) {
    switch (f.kind) {
        case RcfFormulaKind::True:
            return "[1 = 1]";
        case RcfFormulaKind::False:
            return "[1 = 0]";
        case RcfFormulaKind::Atom: {
            std::string op=" = 0";
            if (f.atom.relation==RcfRelation::GreaterEqualZero)
                op=" >= 0";
            else if (f.atom.relation==RcfRelation::GreaterZero)
                op=" > 0";
            return "["+qepcad_poly(f.atom.polynomial)+op+"]";
        }
        case RcfFormulaKind::Not:
            if (f.children.size()!=1)
                throw std::invalid_argument(
                    "prenex matrix not must have one child");
            return "~"+qepcad_matrix(f.children.front());
        case RcfFormulaKind::And:
        case RcfFormulaKind::Or: {
            const char* op=
                f.kind==RcfFormulaKind::And?" /\\ ":" \\/ ";
            std::ostringstream os;
            os << "[";
            for (std::size_t i=0;i<f.children.size();++i) {
                if (i) os << op;
                os << qepcad_matrix(f.children[i]);
            }
            os << "]";
            return os.str();
        }
        case RcfFormulaKind::Exists:
            throw std::invalid_argument(
                "QEPCAD matrix serializer requires prenex normalization");
    }
    throw std::invalid_argument("unknown RCF formula kind");
}

void collect_stats_public(
    const RcfFormula& f,
    std::size_t& max_degree,
    std::size_t& atom_count) {
    if (f.kind==RcfFormulaKind::Atom) {
        ++atom_count;
        max_degree=std::max(
            max_degree,
            rcf_polynomial_degree(f.atom.polynomial));
    }
    for (const auto& c:f.children)
        collect_stats_public(c,max_degree,atom_count);
}

} // namespace

std::size_t rcf_polynomial_degree(
    const RcfPolynomialExpr& p) {
    switch (p.kind) {
        case RcfPolynomialKind::Rational:
            return 0;
        case RcfPolynomialKind::Variable:
            return 1;
        case RcfPolynomialKind::Add: {
            std::size_t d=0;
            for (const auto& x:p.args)
                d=std::max(d,rcf_polynomial_degree(x));
            return d;
        }
        case RcfPolynomialKind::Multiply: {
            std::size_t d=0;
            for (const auto& x:p.args)
                d+=rcf_polynomial_degree(x);
            return d;
        }
        case RcfPolynomialKind::Negate:
            if (p.args.size()!=1)
                throw std::invalid_argument(
                    "RCF polynomial negation must have one child");
            return rcf_polynomial_degree(p.args.front());
    }
    throw std::invalid_argument("unknown RCF polynomial kind");
}

QepcadRcfProgram lower_to_qepcad_rcf_program(
    const ExactSemialgebraicProgram& source) {
    const auto validation=verify_exact_semialgebraic_program(source);
    if (validation.status!=ExactSemialgebraicProgramValidationStatus::Valid)
        throw std::invalid_argument(
            "cannot build QEPCAD RCF program from invalid Phase 6B2B source");
    return RcfLowerer(source).run();
}

QepcadRcfProgram compile_qepcad_rcf_program(
    const NetworkClosureInput& input) {
    return lower_to_qepcad_rcf_program(
        compile_exact_semialgebraic_program(input));
}

std::string serialize_rcf_polynomial_expr(
    const RcfPolynomialExpr& p) {
    std::ostringstream os;
    os << "{\"kind\":\"" << poly_kind_text(p.kind) << "\"";
    if (p.kind==RcfPolynomialKind::Rational)
        os << ",\"value\":\"" << rational_text(p.rational) << "\"";
    if (p.kind==RcfPolynomialKind::Variable)
        os << ",\"variable\":" << p.variable;
    if (!p.args.empty()) {
        os << ",\"args\":[";
        for (std::size_t i=0;i<p.args.size();++i) {
            if (i) os << ",";
            os << serialize_rcf_polynomial_expr(p.args[i]);
        }
        os << "]";
    }
    os << "}";
    return os.str();
}

std::string serialize_rcf_formula(
    const RcfFormula& f) {
    std::ostringstream os;
    os << "{\"kind\":\"" << formula_kind_text(f.kind) << "\"";
    if (f.kind==RcfFormulaKind::Atom) {
        os << ",\"relation\":\"" << relation_text(f.atom.relation)
           << "\",\"polynomial\":"
           << serialize_rcf_polynomial_expr(f.atom.polynomial);
    }
    if (f.kind==RcfFormulaKind::Exists) {
        os << ",\"variables\":";
        append_size_vector(os,f.quantified_variables);
    }
    if (!f.children.empty()) {
        os << ",\"children\":[";
        for (std::size_t i=0;i<f.children.size();++i) {
            if (i) os << ",";
            os << serialize_rcf_formula(f.children[i]);
        }
        os << "]";
    }
    os << "}";
    return os.str();
}

std::string serialize_qepcad_rcf_program(
    const QepcadRcfProgram& p) {
    std::ostringstream os;
    os << "{\"schema\":\"" << p.schema
       << "\",\"schema_version\":\"" << p.schema_version
       << "\",\"backend_contract\":\"" << p.backend_contract
       << "\",\"theorem_id\":\"" << p.theorem_id
       << "\",\"source_digest\":{\"algorithm\":\""
       << p.source_digest.algorithm
       << "\",\"value\":\"" << p.source_digest.value
       << "\"},\"trig_denominator_limit\":" << kMaximumTrigDenominator
       << ",\"source_program_digest\":\"" << p.source_program_digest
       << "\",\"input\":";
    append_input(os,p.input);

    os << ",\"variables\":[";
    for (std::size_t i=0;i<p.variables.size();++i) {
        if (i) os << ",";
        const auto& v=p.variables[i];
        os << "{\"index\":" << v.index
           << ",\"name\":\"" << json_escape(v.name)
           << "\",\"role\":\"" << variable_role_text(v.role) << "\"";
        if (v.role==QeVariableRole::AlgebraicCosine ||
            v.role==QeVariableRole::AlgebraicSine)
            os << ",\"phase_pi\":\"" << rational_string(v.phase_pi) << "\"";
        if (v.role==QeVariableRole::SourceOuter ||
            v.role==QeVariableRole::SourcePairLocal)
            os << ",\"source_variable\":" << v.source_variable;
        if (v.role==QeVariableRole::SourcePairLocal)
            os << ",\"pair\":[" << v.pair_e << "," << v.pair_f << "]";
        os << "}";
    }
    os << "],\"coefficient_variables\":";
    append_size_vector(os,p.coefficient_variables);

    os << ",\"trig_definitions\":[";
    for (std::size_t i=0;i<p.trig_definitions.size();++i) {
        if (i) os << ",";
        const auto& t=p.trig_definitions[i];
        os << "{\"phase_pi\":\"" << rational_string(t.phase_pi)
           << "\",\"cosine_variable\":" << t.cosine_variable
           << ",\"sine_variable\":" << t.sine_variable
           << ",\"cosine_polynomial\":[";
        for (std::size_t j=0;j<t.cosine_polynomial.coefficients.size();++j) {
            if (j) os << ",";
            os << "\"" << t.cosine_polynomial.coefficients[j] << "\"";
        }
        os << "],\"cosine_lower\":\"" << rational_text(t.cosine_lower)
           << "\",\"cosine_upper\":\"" << rational_text(t.cosine_upper)
           << "\",\"sine_sign\":" << t.sine_sign << "}";
    }
    os << "],\"operation_definitions\":[";
    for (std::size_t i=0;i<p.operation_definitions.size();++i) {
        if (i) os << ",";
        const auto& d=p.operation_definitions[i];
        os << "{\"variable\":" << d.variable
           << ",\"source_expression\":\""
           << json_escape(d.source_expression) << "\"}";
    }
    os << "],\"maximum_total_degree\":" << p.maximum_total_degree
       << ",\"rational_atom_count\":" << p.rational_atom_count
       << ",\"formula\":" << serialize_rcf_formula(p.formula)
       << "}";
    return os.str();
}

std::string qepcad_rcf_program_digest(
    const QepcadRcfProgram& program) {
    return sha256_hex(serialize_qepcad_rcf_program(program));
}

std::string export_qepcad_input(
    const QepcadRcfProgram& program) {
    const auto validation=verify_qepcad_rcf_program(program);
    if (validation.status!=QepcadRcfValidationStatus::Valid)
        throw std::invalid_argument(
            "cannot export invalid QEPCAD RCF program");

    auto pn=prenex(program.formula,false);

    std::set<std::size_t> seen;
    for (const auto& q:pn.quantifiers) {
        if (q.variable>=program.variables.size())
            throw std::invalid_argument(
                "QEPCAD prenex variable index out of range");
        if (!seen.insert(q.variable).second)
            throw std::invalid_argument(
                "QEPCAD prenex variable quantified more than once");
    }
    if (seen.size()!=program.variables.size())
        throw std::invalid_argument(
            "QEPCAD decision program does not quantify every variable");

    std::ostringstream os;
    os << "[ PCG exact fixed-turn realizability "
       << qepcad_rcf_program_digest(program) << " ]\n";
    os << "(";
    for (std::size_t i=0;i<pn.quantifiers.size();++i) {
        if (i) os << ",";
        os << "v" << pn.quantifiers[i].variable;
    }
    os << ")\n";
    os << "0\n";
    for (const auto& q:pn.quantifiers)
        os << "(" << (q.existential?"E":"A")
           << " v" << q.variable << ")";
    os << qepcad_matrix(pn.matrix) << ".\n";
    os << "finish\n";
    return os.str();
}

QepcadRcfValidation verify_qepcad_rcf_program(
    const QepcadRcfProgram& program) {
    try {
        if (!metadata_ok(program))
            return {
                QepcadRcfValidationStatus::InvalidMetadata,
                "schema, backend contract, theorem id, or source digest mismatch"};

        const auto source=compile_exact_semialgebraic_program(program.input);
        if (program.source_program_digest!=
            exact_semialgebraic_program_digest(source))
            return {
                QepcadRcfValidationStatus::InvalidSourceBinding,
                "Phase 6B2B source program digest mismatch"};

        for (const auto& t:program.trig_definitions) {
            if (t.cosine_variable>=program.variables.size() ||
                t.sine_variable>=program.variables.size())
                return {
                    QepcadRcfValidationStatus::InvalidAlgebraicIsolation,
                    "trig isolation variable index out of range"};
            if (!(t.cosine_lower<t.cosine_upper))
                return {
                    QepcadRcfValidationStatus::InvalidAlgebraicIsolation,
                    "trig isolation interval is not ordered"};
            if (root_count_open(
                    t.cosine_polynomial,
                    t.cosine_lower,
                    t.cosine_upper)!=1)
                return {
                    QepcadRcfValidationStatus::InvalidAlgebraicIsolation,
                    "cosine interval does not isolate exactly one root"};
            if (t.sine_sign!=sin_pi_sign(t.phase_pi))
                return {
                    QepcadRcfValidationStatus::InvalidAlgebraicIsolation,
                    "sine sign disagrees with exact rational-pi phase"};
        }

        std::size_t max_degree=0,atom_count=0;
        collect_stats_public(
            program.formula,max_degree,atom_count);
        if (max_degree!=program.maximum_total_degree)
            return {
                QepcadRcfValidationStatus::DegreeMismatch,
                "declared rational RCF degree disagrees with formula"};
        if (atom_count!=program.rational_atom_count)
            return {
                QepcadRcfValidationStatus::InvalidStructure,
                "rational atom count mismatch"};

        const auto expected=lower_to_qepcad_rcf_program(source);
        if (serialize_qepcad_rcf_program(program)!=
            serialize_qepcad_rcf_program(expected))
            return {
                QepcadRcfValidationStatus::InvalidStructure,
                "QEPCAD RCF program differs from deterministic recompilation"};

        return {QepcadRcfValidationStatus::Valid,""};
    } catch (const std::exception& e) {
        return {
            QepcadRcfValidationStatus::InvalidStructure,
            e.what()};
    } catch (...) {
        return {
            QepcadRcfValidationStatus::InvalidStructure,
            "unknown QEPCAD RCF verification failure"};
    }
}

} // namespace pcg
