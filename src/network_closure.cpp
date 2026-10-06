#include <pcg/network_closure.hpp>

#include <algorithm>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace pcg {
namespace {

constexpr unsigned kMaxCyclotomicOrder = 256;
const char* kSchema = "pcg-network-closure-certificate";
const char* kSchemaVersion = "1.0";
const char* kTheoremId = "PCG-NETWORK-CLOSURE-THM-001";
const char* kSourceDigest = "6d1f8fc39b3bb34385956c4d2c0e77e222260f60d2f59fde860d6c1964aa8546";

BigInt abs_big(BigInt x) {
    return x < 0 ? -x : x;
}

BigInt gcd_big(BigInt a, BigInt b) {
    a=abs_big(std::move(a));
    b=abs_big(std::move(b));
    while (b != 0) {
        BigInt r=a%b;
        a=std::move(b);
        b=std::move(r);
    }
    return a;
}

std::string rational_text(const Rational& r) {
    std::ostringstream os;
    os << r.numerator() << "/" << r.denominator();
    return os.str();
}

void validate_input(const NetworkClosureInput& input) {
    if (input.vertex_count == 0)
        throw std::invalid_argument("network closure requires at least one quotient vertex");
    if (input.trace_vertices.empty())
        throw std::invalid_argument("network closure requires at least one source interval");
    if (input.trace_vertices.size() != input.turns.size())
        throw std::invalid_argument("trace and turn counts must match");

    std::vector<bool> seen(input.vertex_count,false);
    for (std::size_t v : input.trace_vertices) {
        if (v >= input.vertex_count)
            throw std::invalid_argument("trace vertex index out of range");
        seen[v]=true;
    }
    for (bool x : seen)
        if (!x) throw std::invalid_argument("every declared quotient vertex must occur in the cyclic trace");
}

struct Dsu {
    explicit Dsu(std::size_t n) : p(n), rank(n,0) {
        for (std::size_t i=0;i<n;++i) p[i]=i;
    }
    std::size_t find(std::size_t x) {
        if (p[x]!=x) p[x]=find(p[x]);
        return p[x];
    }
    bool unite(std::size_t a,std::size_t b) {
        a=find(a); b=find(b);
        if (a==b) return false;
        if (rank[a]<rank[b]) std::swap(a,b);
        p[b]=a;
        if (rank[a]==rank[b]) ++rank[a];
        return true;
    }
    std::vector<std::size_t> p;
    std::vector<unsigned> rank;
};

NetworkCycleBasis cycle_basis_impl(const NetworkClosureInput& input) {
    validate_input(input);
    const std::size_t n=input.vertex_count;
    const std::size_t m=input.trace_vertices.size();

    struct Edge { std::size_t tail,head; bool tree{false}; };
    std::vector<Edge> edges(m);
    for (std::size_t i=0;i<m;++i)
        edges[i]={input.trace_vertices[i],input.trace_vertices[(i+1)%m],false};

    Dsu dsu(n);
    std::vector<std::vector<std::pair<std::size_t,std::size_t>>> tree(n);
    std::vector<std::size_t> non_tree;
    for (std::size_t i=0;i<m;++i) {
        const auto& e=edges[i];
        if (e.tail!=e.head && dsu.unite(e.tail,e.head)) {
            edges[i].tree=true;
            tree[e.tail].push_back({e.head,i});
            tree[e.head].push_back({e.tail,i});
        } else {
            non_tree.push_back(i);
        }
    }

    if (non_tree.size() != m-n+1)
        throw std::runtime_error("cyclic trace quotient graph is unexpectedly disconnected");

    NetworkCycleBasis out;
    out.rows.reserve(non_tree.size());

    for (std::size_t eidx : non_tree) {
        std::vector<int> row(m,0);
        row[eidx]=+1;
        const auto& e=edges[eidx];

        if (e.tail != e.head) {
            std::vector<std::size_t> parent_v(n,n);
            std::queue<std::size_t> q;
            parent_v[e.head]=e.head;
            q.push(e.head);
            while(!q.empty() && parent_v[e.tail]==n) {
                const auto v=q.front(); q.pop();
                for (const auto& [w,te] : tree[v]) {
                    if (parent_v[w]!=n) continue;
                    parent_v[w]=v;
                    q.push(w);
                }
            }
            if (parent_v[e.tail]==n)
                throw std::runtime_error("spanning-tree path missing");

            std::size_t cur=e.tail;
            // parent chain runs from tail back toward head. We need the
            // actual cycle traversal head -> ... -> tail, so collect then reverse.
            std::vector<std::pair<std::size_t,std::size_t>> steps;
            while (cur != e.head) {
                const std::size_t pv=parent_v[cur];
                steps.push_back({pv,cur});
                cur=pv;
            }
            std::reverse(steps.begin(),steps.end());

            for (const auto& [a,b] : steps) {
                std::size_t te=m;
                for (const auto& [w,idx] : tree[a]) {
                    if (w==b) { te=idx; break; }
                }
                if (te==m) throw std::runtime_error("tree edge lookup failed");
                row[te]=(edges[te].tail==a && edges[te].head==b) ? +1 : -1;
            }
        }
        out.rows.push_back(std::move(row));
    }

    out.chord_phase_pi.resize(m);
    PiRational phase{0,1};
    for (std::size_t i=0;i<m;++i) {
        out.chord_phase_pi[i]=phase + input.turns[i].pi/2;
        phase=phase+input.turns[i].pi;
    }
    return out;
}

using IntPoly = std::vector<BigInt>;
using RatPoly = std::vector<Rational>;

void trim(IntPoly& p) {
    while (p.size()>1 && p.back()==0) p.pop_back();
}
void trim(RatPoly& p) {
    while (p.size()>1 && p.back()==Rational(BigInt{0})) p.pop_back();
}

IntPoly divide_exact_int(IntPoly a, const IntPoly& b) {
    IntPoly divisor=b;
    trim(a); trim(divisor);
    if (divisor.empty() || divisor.back()==0)
        throw std::runtime_error("zero polynomial divisor");
    if (a.size()<divisor.size()) throw std::runtime_error("non-exact cyclotomic division");

    IntPoly q(a.size()-divisor.size()+1,BigInt{0});
    while (a.size()>=divisor.size() && !(a.size()==1 && a[0]==0)) {
        const std::size_t shift=a.size()-divisor.size();
        if (divisor.back()!=1)
            throw std::runtime_error("cyclotomic divisor must be monic");
        const BigInt coeff=a.back();
        q[shift]=coeff;
        for (std::size_t j=0;j<divisor.size();++j)
            a[shift+j]-=coeff*divisor[j];
        trim(a);
    }
    if (!(a.size()==1 && a[0]==0))
        throw std::runtime_error("cyclotomic polynomial division left remainder");
    trim(q);
    return q;
}

IntPoly cyclotomic(unsigned n, std::map<unsigned,IntPoly>& memo) {
    auto it=memo.find(n);
    if (it!=memo.end()) return it->second;

    IntPoly p(n+1,BigInt{0});
    p[0]=-1;
    p[n]=1;
    for (unsigned d=1;d<n;++d) {
        if (n%d==0)
            p=divide_exact_int(std::move(p),cyclotomic(d,memo));
    }
    trim(p);
    memo[n]=p;
    return p;
}

struct Alg {
    std::vector<Rational> c;
};

struct CyclotomicField {
    unsigned order{};
    IntPoly phi;
    std::size_t degree{};

    explicit CyclotomicField(unsigned n) : order(n) {
        std::map<unsigned,IntPoly> memo;
        phi=cyclotomic(n,memo);
        degree=phi.size()-1;
        if (degree==0) throw std::runtime_error("invalid cyclotomic degree");
    }

    Alg zero() const { return Alg{std::vector<Rational>(degree,Rational(BigInt{0}))}; }
    Alg one() const {
        Alg a=zero();
        a.c[0]=Rational(BigInt{1});
        return a;
    }

    bool is_zero(const Alg& a) const {
        for (const auto& x : a.c)
            if (x != Rational(BigInt{0})) return false;
        return true;
    }

    bool equal(const Alg& a,const Alg& b) const {
        return a.c==b.c;
    }

    Alg reduce(RatPoly p) const {
        if (p.empty()) p.push_back(Rational(BigInt{0}));
        while (p.size()>degree) {
            const std::size_t idx=p.size()-1;
            const Rational factor=p.back();
            if (factor != Rational(BigInt{0})) {
                const std::size_t shift=idx-degree;
                for (std::size_t j=0;j<phi.size();++j)
                    p[shift+j]-=factor*Rational(phi[j]);
            }
            trim(p);
        }
        p.resize(degree,Rational(BigInt{0}));
        return Alg{std::move(p)};
    }

    Alg monomial(unsigned exponent) const {
        exponent%=order;
        RatPoly p(exponent+1,Rational(BigInt{0}));
        p[exponent]=Rational(BigInt{1});
        return reduce(std::move(p));
    }

    Alg add(const Alg& a,const Alg& b) const {
        Alg r=zero();
        for (std::size_t i=0;i<degree;++i) r.c[i]=a.c[i]+b.c[i];
        return r;
    }
    Alg sub(const Alg& a,const Alg& b) const {
        Alg r=zero();
        for (std::size_t i=0;i<degree;++i) r.c[i]=a.c[i]-b.c[i];
        return r;
    }
    Alg neg(const Alg& a) const {
        Alg r=zero();
        for (std::size_t i=0;i<degree;++i) r.c[i]=-a.c[i];
        return r;
    }
    Alg scale(const Alg& a,const Rational& s) const {
        Alg r=zero();
        for (std::size_t i=0;i<degree;++i) r.c[i]=a.c[i]*s;
        return r;
    }
    Alg mul(const Alg& a,const Alg& b) const {
        RatPoly p(2*degree-1,Rational(BigInt{0}));
        for (std::size_t i=0;i<degree;++i)
            for (std::size_t j=0;j<degree;++j)
                p[i+j]+=a.c[i]*b.c[j];
        return reduce(std::move(p));
    }

    static std::pair<RatPoly,RatPoly> divrem_poly(RatPoly a, RatPoly b) {
        trim(a); trim(b);
        if (b.size()==1 && b[0]==Rational(BigInt{0}))
            throw std::runtime_error("zero polynomial divisor");
        RatPoly q(a.size()>=b.size()?a.size()-b.size()+1:1,Rational(BigInt{0}));
        while (!(a.size()==1 && a[0]==Rational(BigInt{0})) && a.size()>=b.size()) {
            const std::size_t shift=a.size()-b.size();
            const Rational coeff=a.back()/b.back();
            q[shift]+=coeff;
            for (std::size_t j=0;j<b.size();++j)
                a[shift+j]-=coeff*b[j];
            trim(a);
        }
        trim(q);
        return {q,a};
    }

    static RatPoly poly_sub(const RatPoly& a,const RatPoly& b) {
        RatPoly r(std::max(a.size(),b.size()),Rational(BigInt{0}));
        for (std::size_t i=0;i<a.size();++i) r[i]+=a[i];
        for (std::size_t i=0;i<b.size();++i) r[i]-=b[i];
        trim(r);
        return r;
    }

    static RatPoly poly_mul(const RatPoly& a,const RatPoly& b) {
        RatPoly r(a.size()+b.size()-1,Rational(BigInt{0}));
        for (std::size_t i=0;i<a.size();++i)
            for (std::size_t j=0;j<b.size();++j)
                r[i+j]+=a[i]*b[j];
        trim(r);
        return r;
    }

    Alg inverse(const Alg& a) const {
        if (is_zero(a)) throw std::runtime_error("division by zero algebraic element");

        RatPoly r0(phi.size(),Rational(BigInt{0}));
        for (std::size_t i=0;i<phi.size();++i) r0[i]=Rational(phi[i]);
        RatPoly r1=a.c; trim(r1);
        RatPoly t0{Rational(BigInt{0})};
        RatPoly t1{Rational(BigInt{1})};

        while (!(r1.size()==1 && r1[0]==Rational(BigInt{0}))) {
            auto [q,r]=divrem_poly(r0,r1);
            RatPoly next_t=poly_sub(t0,poly_mul(q,t1));
            r0=std::move(r1);
            r1=std::move(r);
            t0=std::move(t1);
            t1=std::move(next_t);
        }

        trim(r0);
        if (r0.size()!=1 || r0[0]==Rational(BigInt{0}))
            throw std::runtime_error("algebraic inverse gcd failure");
        for (auto& x : t0) x/=r0[0];
        return reduce(std::move(t0));
    }

    Alg divide(const Alg& a,const Alg& b) const {
        return mul(a,inverse(b));
    }

    Alg conjugate(const Alg& a) const {
        Alg r=zero();
        for (std::size_t j=0;j<degree;++j) {
            if (a.c[j]==Rational(BigInt{0})) continue;
            const unsigned exp=(j==0)?0U:(order-static_cast<unsigned>(j%order))%order;
            r=add(r,scale(monomial(exp),a.c[j]));
        }
        return r;
    }

    bool is_real(const Alg& a) const {
        return equal(a,conjugate(a));
    }
};

detail::Interval rational_interval_local(const Rational& q) {
    using Dec100 = boost::multiprecision::number<boost::multiprecision::cpp_dec_float<100>>;
    const Dec100 approx=Dec100(q.numerator())/Dec100(q.denominator());
    const double d=approx.convert_to<double>();
    return detail::Interval(
        std::nextafter(d,-std::numeric_limits<double>::infinity()),
        std::nextafter(d, std::numeric_limits<double>::infinity()));
}

enum class SignCert { Negative, Zero, Positive, Indeterminate };

SignCert certified_real_sign(const CyclotomicField& field,const Alg& a) {
    if (field.is_zero(a)) return SignCert::Zero;
    if (!field.is_real(a)) return SignCert::Indeterminate;

    detail::Interval acc(0.0);
    for (std::size_t j=0;j<field.degree;++j) {
        if (a.c[j]==Rational(BigInt{0})) continue;
        const PiRational phase{BigInt{2}*BigInt{j},BigInt{field.order}};
        acc += rational_interval_local(a.c[j]) * certified_cos_pi(phase).interval();
    }
    if (acc.lower()>0.0) return SignCert::Positive;
    if (acc.upper()<0.0) return SignCert::Negative;
    return SignCert::Indeterminate;
}

struct ExactSystem {
    bool supported{false};
    unsigned order{};
    NetworkCycleBasis topology;
    CyclotomicField* field_ptr{nullptr};
    std::vector<std::vector<Alg>> matrix;
};

struct OwnedExactSystem {
    bool supported{false};
    unsigned order{};
    NetworkCycleBasis topology;
    CyclotomicField field{1};
    std::vector<std::vector<Alg>> matrix;

    OwnedExactSystem() : field(1) {}
    explicit OwnedExactSystem(unsigned n) : supported(true), order(n), field(n) {}
};

OwnedExactSystem build_exact_system(const NetworkClosureInput& input) {
    NetworkCycleBasis topology=cycle_basis_impl(input);

    BigInt L=1;
    for (const auto& a : topology.chord_phase_pi) {
        const BigInt d=a.value.denominator();
        const BigInt g=gcd_big(L,d);
        const BigInt candidate=(L/g)*d;
        if (candidate*4 > kMaxCyclotomicOrder) {
            OwnedExactSystem out;
            out.topology=std::move(topology);
            return out;
        }
        L=candidate;
    }

    const unsigned Lu=L.convert_to<unsigned>();
    const unsigned order=4*Lu;
    OwnedExactSystem out(order);
    out.topology=std::move(topology);

    const std::size_t beta=out.topology.rows.size();
    const std::size_t m=input.turns.size();
    out.matrix.assign(2*beta,std::vector<Alg>(m,out.field.zero()));

    const Alg imag_unit=out.field.monomial(Lu);
    const Rational half(BigInt{1},BigInt{2});

    for (std::size_t e=0;e<m;++e) {
        const auto& alpha=out.topology.chord_phase_pi[e].value;
        const BigInt scaled=alpha.numerator()*(L/alpha.denominator());
        BigInt exp_big=BigInt{2}*scaled;
        exp_big%=order;
        if (exp_big<0) exp_big+=order;
        const unsigned exp=exp_big.convert_to<unsigned>();
        const unsigned negexp=(order-exp)%order;

        const Alg z=out.field.monomial(exp);
        const Alg zbar=out.field.monomial(negexp);
        const Alg cosine=out.field.scale(out.field.add(z,zbar),half);
        const Alg sine=out.field.scale(
            out.field.mul(out.field.sub(z,zbar),out.field.neg(imag_unit)),
            half);

        if (!out.field.is_real(cosine) || !out.field.is_real(sine))
            throw std::runtime_error("trigonometric cyclotomic construction is not exactly real");

        for (std::size_t k=0;k<beta;++k) {
            const int coeff=out.topology.rows[k][e];
            if (coeff==0) continue;
            const Rational s(BigInt{coeff});
            out.matrix[k][e]=out.field.scale(cosine,s);
            out.matrix[beta+k][e]=out.field.scale(sine,s);
        }
    }
    return out;
}

struct RrefResult {
    std::vector<std::vector<Alg>> r;
    std::vector<std::size_t> pivot_cols;
};

RrefResult rref(const CyclotomicField& field,
                std::vector<std::vector<Alg>> a) {
    if (a.empty()) return {std::move(a),{}};
    const std::size_t rows=a.size();
    const std::size_t cols=a[0].size();
    std::vector<std::size_t> pivots;
    std::size_t row=0;

    for (std::size_t col=0;col<cols && row<rows;++col) {
        std::size_t pivot=row;
        while (pivot<rows && field.is_zero(a[pivot][col])) ++pivot;
        if (pivot==rows) continue;
        if (pivot!=row) std::swap(a[pivot],a[row]);

        const Alg inv=field.inverse(a[row][col]);
        for (std::size_t j=0;j<cols;++j)
            a[row][j]=field.mul(a[row][j],inv);

        for (std::size_t i=0;i<rows;++i) {
            if (i==row || field.is_zero(a[i][col])) continue;
            const Alg factor=a[i][col];
            for (std::size_t j=0;j<cols;++j)
                a[i][j]=field.sub(a[i][j],field.mul(factor,a[row][j]));
        }
        pivots.push_back(col);
        ++row;
    }
    return {std::move(a),std::move(pivots)};
}

std::vector<Alg> rigid_kernel_from_rref(const CyclotomicField& field,
                                        const RrefResult& rr,
                                        std::size_t cols) {
    std::vector<bool> pivot(cols,false);
    for (auto c : rr.pivot_cols) pivot[c]=true;
    std::size_t free_col=cols;
    for (std::size_t c=0;c<cols;++c)
        if (!pivot[c]) { free_col=c; break; }
    if (free_col==cols) throw std::runtime_error("no free column in rigid kernel");

    std::vector<Alg> x(cols,field.zero());
    x[free_col]=field.one();
    for (std::size_t r=0;r<rr.pivot_cols.size();++r) {
        const std::size_t pc=rr.pivot_cols[r];
        x[pc]=field.neg(rr.r[r][free_col]);
    }
    return x;
}

bool matrix_times_vector_zero(const CyclotomicField& field,
                              const std::vector<std::vector<Alg>>& a,
                              const std::vector<Alg>& x) {
    for (const auto& row : a) {
        Alg sum=field.zero();
        for (std::size_t j=0;j<x.size();++j)
            sum=field.add(sum,field.mul(row[j],x[j]));
        if (!field.is_zero(sum)) return false;
    }
    return true;
}

std::vector<AlgebraicChordValue> encode_vector(const std::vector<Alg>& x) {
    std::vector<AlgebraicChordValue> out;
    out.reserve(x.size());
    for (const auto& a : x) out.push_back(AlgebraicChordValue{a.c});
    return out;
}

std::vector<Alg> decode_vector(const CyclotomicField& field,
                               const std::vector<AlgebraicChordValue>& x) {
    std::vector<Alg> out;
    out.reserve(x.size());
    for (const auto& v : x) {
        if (v.coefficients.size()!=field.degree)
            throw std::invalid_argument("noncanonical algebraic certificate coefficient count");
        out.push_back(Alg{v.coefficients});
    }
    return out;
}

struct SignSummary {
    bool indeterminate{false};
    bool has_zero{false};
    bool has_pos{false};
    bool has_neg{false};
};

SignSummary summarize_signs(const CyclotomicField& field,const std::vector<Alg>& x) {
    SignSummary s;
    for (const auto& a : x) {
        const auto sign=certified_real_sign(field,a);
        if (sign==SignCert::Indeterminate) s.indeterminate=true;
        if (sign==SignCert::Zero) s.has_zero=true;
        if (sign==SignCert::Positive) s.has_pos=true;
        if (sign==SignCert::Negative) s.has_neg=true;
    }
    return s;
}


std::vector<std::vector<Alg>> kernel_basis_from_rref(const CyclotomicField& field,
                                                     const RrefResult& rr,
                                                     std::size_t cols) {
    std::vector<bool> pivot(cols,false);
    for (auto c : rr.pivot_cols) pivot[c]=true;

    std::vector<std::size_t> free_cols;
    for (std::size_t c=0;c<cols;++c)
        if (!pivot[c]) free_cols.push_back(c);

    std::vector<std::vector<Alg>> basis;
    basis.reserve(free_cols.size());
    for (const std::size_t free_col : free_cols) {
        std::vector<Alg> x(cols,field.zero());
        x[free_col]=field.one();
        for (std::size_t r=0;r<rr.pivot_cols.size();++r) {
            const std::size_t pc=rr.pivot_cols[r];
            x[pc]=field.neg(rr.r[r][free_col]);
        }
        basis.push_back(std::move(x));
    }
    return basis;
}

struct LinearInequality {
    std::vector<Alg> coeff;
    Alg rhs;
};

enum class ConeState {
    Feasible,
    Infeasible,
    Indeterminate,
    ComplexityLimit
};

struct ConeResult {
    ConeState state{ConeState::Indeterminate};
    std::vector<Alg> parameters;
};

constexpr std::size_t kMaxEliminationInequalities = 8192;

SignCert compare_real(const CyclotomicField& field,const Alg& a,const Alg& b) {
    return certified_real_sign(field,field.sub(a,b));
}

Alg eval_linear(const CyclotomicField& field,
                const std::vector<Alg>& coeff,
                const std::vector<Alg>& x) {
    Alg sum=field.zero();
    for (std::size_t i=0;i<coeff.size();++i)
        sum=field.add(sum,field.mul(coeff[i],x[i]));
    return sum;
}

ConeResult solve_linear_inequalities(const CyclotomicField& field,
                                     const std::vector<LinearInequality>& ineq,
                                     std::size_t vars) {
    if (vars==0) {
        for (const auto& q : ineq) {
            if (!q.coeff.empty()) return {ConeState::Indeterminate,{}};
            const auto sign=certified_real_sign(field,q.rhs);
            if (sign==SignCert::Positive) return {ConeState::Infeasible,{}};
            if (sign==SignCert::Indeterminate) return {ConeState::Indeterminate,{}};
        }
        return {ConeState::Feasible,{}};
    }

    const std::size_t xcol=vars-1;
    std::vector<const LinearInequality*> pos,neg,zero;
    pos.reserve(ineq.size());
    neg.reserve(ineq.size());
    zero.reserve(ineq.size());

    for (const auto& q : ineq) {
        if (q.coeff.size()!=vars) return {ConeState::Indeterminate,{}};
        const auto s=certified_real_sign(field,q.coeff[xcol]);
        if (s==SignCert::Positive) pos.push_back(&q);
        else if (s==SignCert::Negative) neg.push_back(&q);
        else if (s==SignCert::Zero) zero.push_back(&q);
        else return {ConeState::Indeterminate,{}};
    }

    if (zero.size() > kMaxEliminationInequalities)
        return {ConeState::ComplexityLimit,{}};
    if (!pos.empty() && !neg.empty() &&
        pos.size() > kMaxEliminationInequalities / neg.size())
        return {ConeState::ComplexityLimit,{}};

    std::vector<LinearInequality> reduced;
    reduced.reserve(zero.size()+pos.size()*neg.size());

    for (const auto* q : zero) {
        LinearInequality r;
        r.coeff.assign(q->coeff.begin(),q->coeff.begin()+static_cast<std::ptrdiff_t>(xcol));
        r.rhs=q->rhs;
        reduced.push_back(std::move(r));
    }

    for (const auto* p : pos) {
        const Alg minus_n_dummy=field.zero();
        (void)minus_n_dummy;
        for (const auto* n : neg) {
            const Alg mp=field.neg(n->coeff[xcol]); // -negative coefficient > 0
            const Alg pp=p->coeff[xcol];            // positive coefficient
            LinearInequality r;
            r.coeff.resize(xcol,field.zero());
            for (std::size_t j=0;j<xcol;++j) {
                r.coeff[j]=field.add(
                    field.mul(mp,p->coeff[j]),
                    field.mul(pp,n->coeff[j]));
            }
            r.rhs=field.add(field.mul(mp,p->rhs),field.mul(pp,n->rhs));
            reduced.push_back(std::move(r));
            if (reduced.size()>kMaxEliminationInequalities)
                return {ConeState::ComplexityLimit,{}};
        }
    }

    auto sub=solve_linear_inequalities(field,reduced,vars-1);
    if (sub.state!=ConeState::Feasible) return sub;

    bool have_lower=false,have_upper=false;
    Alg lower=field.zero(),upper=field.zero();

    for (const auto& q : ineq) {
        const auto ax_sign=certified_real_sign(field,q.coeff[xcol]);
        if (ax_sign==SignCert::Zero) {
            std::vector<Alg> prefix(q.coeff.begin(),q.coeff.begin()+static_cast<std::ptrdiff_t>(xcol));
            const Alg lhs=eval_linear(field,prefix,sub.parameters);
            const auto ok=compare_real(field,lhs,q.rhs);
            if (ok==SignCert::Negative) return {ConeState::Infeasible,{}};
            if (ok==SignCert::Indeterminate) return {ConeState::Indeterminate,{}};
            continue;
        }
        if (ax_sign==SignCert::Indeterminate)
            return {ConeState::Indeterminate,{}};

        std::vector<Alg> prefix(q.coeff.begin(),q.coeff.begin()+static_cast<std::ptrdiff_t>(xcol));
        const Alg rest=eval_linear(field,prefix,sub.parameters);
        const Alg bound=field.divide(field.sub(q.rhs,rest),q.coeff[xcol]);

        if (ax_sign==SignCert::Positive) {
            if (!have_lower) {
                lower=bound;
                have_lower=true;
            } else {
                const auto cmp=compare_real(field,bound,lower);
                if (cmp==SignCert::Positive) lower=bound;
                else if (cmp==SignCert::Indeterminate) return {ConeState::Indeterminate,{}};
            }
        } else {
            if (!have_upper) {
                upper=bound;
                have_upper=true;
            } else {
                const auto cmp=compare_real(field,bound,upper);
                if (cmp==SignCert::Negative) upper=bound;
                else if (cmp==SignCert::Indeterminate) return {ConeState::Indeterminate,{}};
            }
        }
    }

    Alg x=field.zero();
    if (have_lower && have_upper) {
        const auto gap=compare_real(field,upper,lower);
        if (gap==SignCert::Negative) return {ConeState::Infeasible,{}};
        if (gap==SignCert::Indeterminate) return {ConeState::Indeterminate,{}};
        x=field.scale(field.add(lower,upper),Rational(BigInt{1},BigInt{2}));
    } else if (have_lower) {
        x=lower;
    } else if (have_upper) {
        x=upper;
    }

    std::vector<Alg> witness=sub.parameters;
    witness.push_back(x);

    for (const auto& q : ineq) {
        const Alg lhs=eval_linear(field,q.coeff,witness);
        const auto ok=compare_real(field,lhs,q.rhs);
        if (ok==SignCert::Negative) return {ConeState::Infeasible,{}};
        if (ok==SignCert::Indeterminate) return {ConeState::Indeterminate,{}};
    }
    return {ConeState::Feasible,std::move(witness)};
}

ConeResult solve_general_positive_kernel(const CyclotomicField& field,
                                         const RrefResult& rr,
                                         const std::vector<std::vector<Alg>>& matrix,
                                         std::size_t cols,
                                         std::vector<Alg>* positive_kernel) {
    const auto basis=kernel_basis_from_rref(field,rr,cols);
    if (basis.empty()) return {ConeState::Infeasible,{}};

    std::vector<LinearInequality> inequalities;
    inequalities.reserve(cols);
    for (std::size_t i=0;i<cols;++i) {
        LinearInequality q;
        q.coeff.resize(basis.size(),field.zero());
        for (std::size_t j=0;j<basis.size();++j)
            q.coeff[j]=basis[j][i];
        q.rhs=field.one(); // homogeneity lets strict c_i>0 scale to c_i>=1
        inequalities.push_back(std::move(q));
    }

    auto solved=solve_linear_inequalities(field,inequalities,basis.size());
    if (solved.state!=ConeState::Feasible) return solved;

    std::vector<Alg> kernel(cols,field.zero());
    for (std::size_t j=0;j<basis.size();++j)
        for (std::size_t i=0;i<cols;++i)
            kernel[i]=field.add(kernel[i],field.mul(basis[j][i],solved.parameters[j]));

    if (!matrix_times_vector_zero(field,matrix,kernel))
        return {ConeState::Indeterminate,{}};
    const auto signs=summarize_signs(field,kernel);
    if (signs.indeterminate) return {ConeState::Indeterminate,{}};
    if (signs.has_zero || signs.has_neg || !signs.has_pos)
        return {ConeState::Indeterminate,{}};

    if (positive_kernel) *positive_kernel=std::move(kernel);
    return solved;
}

NetworkClosureCertificate base_certificate(const NetworkClosureInput& input) {
    NetworkClosureCertificate cert(input);
    cert.canonical_input_digest=sha256_hex(canonicalize_network_closure_input(input));
    return cert;
}

NetworkClosureResult result_from_certificate(const NetworkClosureCertificate& cert) {
    NetworkClosureResult out(cert);
    out.status=cert.status;
    out.proof=cert.proof;
    if (cert.proof==NetworkClosureProof::ExactFullRankObstruction)
        out.assurance=ArithmeticAssurance::Exact;
    else if (cert.proof==NetworkClosureProof::CertifiedPositiveRigidKernel ||
             cert.proof==NetworkClosureProof::CertifiedRigidSignObstruction ||
             cert.proof==NetworkClosureProof::CertifiedGeneralPositiveKernel ||
             cert.proof==NetworkClosureProof::CertifiedGeneralConeObstruction)
        out.assurance=ArithmeticAssurance::CertifiedNumerical;
    else
        out.assurance=ArithmeticAssurance::None;
    out.termination=cert.status==NetworkClosureStatus::Indeterminate
        ? (cert.proof==NetworkClosureProof::CyclotomicOrderLimit ||
           cert.proof==NetworkClosureProof::HigherDimensionalKernel ||
           cert.proof==NetworkClosureProof::SignCertificationLimit ||
           cert.proof==NetworkClosureProof::EliminationComplexityLimit
              ? TerminationReason::PrecisionLimit
              : TerminationReason::BackendFailure)
        : TerminationReason::Completed;
    return out;
}

NetworkClosureResult invalid_certificate(const NetworkClosureCertificate& cert) {
    NetworkClosureResult out(cert);
    out.status=NetworkClosureStatus::Indeterminate;
    out.proof=NetworkClosureProof::BackendFailure;
    out.assurance=ArithmeticAssurance::None;
    out.termination=TerminationReason::BackendFailure;
    return out;
}

bool metadata_ok(const NetworkClosureCertificate& c) {
    return c.schema==kSchema &&
           c.schema_version==kSchemaVersion &&
           c.theorem_id==kTheoremId &&
           c.source_digest.algorithm=="SHA-256" &&
           c.source_digest.value==kSourceDigest &&
           c.canonical_input_digest==sha256_hex(canonicalize_network_closure_input(c.input));
}

} // namespace

NetworkCycleBasis build_network_cycle_basis(const NetworkClosureInput& input) {
    return cycle_basis_impl(input);
}

std::string canonicalize_network_closure_input(const NetworkClosureInput& input) {
    validate_input(input);
    std::ostringstream os;
    os << "NETWORK-CLOSURE|vertices=" << input.vertex_count << "|trace=";
    for (std::size_t i=0;i<input.trace_vertices.size();++i) {
        if (i) os << ",";
        os << input.trace_vertices[i];
    }
    os << "|turns=";
    for (std::size_t i=0;i<input.turns.size();++i) {
        if (i) os << ",";
        os << rational_string(input.turns[i].pi);
    }
    return os.str();
}

std::string serialize_network_closure_certificate(const NetworkClosureCertificate& c) {
    std::ostringstream os;
    os << "{\"schema\":\"" << c.schema
       << "\",\"schema_version\":\"" << c.schema_version
       << "\",\"theorem_id\":\"" << c.theorem_id
       << "\",\"source_digest\":\"" << c.source_digest.value
       << "\",\"canonical_input_digest\":\"" << c.canonical_input_digest
       << "\",\"status\":" << static_cast<int>(c.status)
       << ",\"proof\":" << static_cast<int>(c.proof)
       << ",\"cyclotomic_order\":" << c.cyclotomic_order
       << ",\"exact_rank\":" << c.exact_rank
       << ",\"kernel\":[";
    for (std::size_t i=0;i<c.rigid_kernel.size();++i) {
        if (i) os << ",";
        os << "[";
        for (std::size_t j=0;j<c.rigid_kernel[i].coefficients.size();++j) {
            if (j) os << ",";
            os << "\"" << rational_text(c.rigid_kernel[i].coefficients[j]) << "\"";
        }
        os << "]";
    }
    os << "]}";
    return os.str();
}

NetworkClosureResult solve_network_closure(const NetworkClosureInput& input) {
    validate_input(input);
    auto cert=base_certificate(input);

    try {
        auto sys=build_exact_system(input);
        if (!sys.supported) {
            cert.status=NetworkClosureStatus::Indeterminate;
            cert.proof=NetworkClosureProof::CyclotomicOrderLimit;
            return result_from_certificate(cert);
        }

        cert.cyclotomic_order=sys.order;
        const auto rr=rref(sys.field,sys.matrix);
        cert.exact_rank=rr.pivot_cols.size();
        const std::size_t m=input.turns.size();

        if (cert.exact_rank==m) {
            cert.status=NetworkClosureStatus::NotClosed;
            cert.proof=NetworkClosureProof::ExactFullRankObstruction;
        } else if (m-cert.exact_rank==1) {
            auto kernel=rigid_kernel_from_rref(sys.field,rr,m);
            if (!matrix_times_vector_zero(sys.field,sys.matrix,kernel))
                return invalid_certificate(cert);

            const auto signs=summarize_signs(sys.field,kernel);
            if (signs.indeterminate) {
                cert.status=NetworkClosureStatus::Indeterminate;
                cert.proof=NetworkClosureProof::SignCertificationLimit;
                cert.rigid_kernel=encode_vector(kernel);
                return result_from_certificate(cert);
            }

            if (!signs.has_zero && !(signs.has_pos && signs.has_neg)) {
                if (signs.has_neg) {
                    for (auto& x : kernel) x=sys.field.neg(x);
                }
                cert.status=NetworkClosureStatus::Closed;
                cert.proof=NetworkClosureProof::CertifiedPositiveRigidKernel;
                cert.rigid_kernel=encode_vector(kernel);
            } else {
                cert.status=NetworkClosureStatus::NotClosed;
                cert.proof=NetworkClosureProof::CertifiedRigidSignObstruction;
                cert.rigid_kernel=encode_vector(kernel);
            }
        } else {
            std::vector<Alg> kernel;
            const auto cone=solve_general_positive_kernel(
                sys.field,rr,sys.matrix,m,&kernel);

            if (cone.state==ConeState::Feasible) {
                cert.status=NetworkClosureStatus::Closed;
                cert.proof=NetworkClosureProof::CertifiedGeneralPositiveKernel;
                cert.rigid_kernel=encode_vector(kernel);
            } else if (cone.state==ConeState::Infeasible) {
                cert.status=NetworkClosureStatus::NotClosed;
                cert.proof=NetworkClosureProof::CertifiedGeneralConeObstruction;
            } else if (cone.state==ConeState::ComplexityLimit) {
                cert.status=NetworkClosureStatus::Indeterminate;
                cert.proof=NetworkClosureProof::EliminationComplexityLimit;
            } else {
                cert.status=NetworkClosureStatus::Indeterminate;
                cert.proof=NetworkClosureProof::SignCertificationLimit;
            }
        }

        auto checked=verify_network_closure_certificate(cert);
        if (checked.status!=cert.status || checked.proof!=cert.proof)
            return invalid_certificate(cert);
        return result_from_certificate(cert);
    } catch (...) {
        return invalid_certificate(cert);
    }
}

NetworkClosureResult verify_network_closure_certificate(const NetworkClosureCertificate& cert) {
    try {
        validate_input(cert.input);
        if (!metadata_ok(cert)) return invalid_certificate(cert);

        auto sys=build_exact_system(cert.input);
        if (!sys.supported) {
            if (cert.status==NetworkClosureStatus::Indeterminate &&
                cert.proof==NetworkClosureProof::CyclotomicOrderLimit &&
                cert.cyclotomic_order==0 &&
                cert.exact_rank==0 &&
                cert.rigid_kernel.empty())
                return result_from_certificate(cert);
            return invalid_certificate(cert);
        }

        if (cert.cyclotomic_order!=sys.order) return invalid_certificate(cert);
        const auto rr=rref(sys.field,sys.matrix);
        const std::size_t rank=rr.pivot_cols.size();
        const std::size_t m=cert.input.turns.size();
        if (cert.exact_rank!=rank) return invalid_certificate(cert);

        if (cert.status==NetworkClosureStatus::NotClosed &&
            cert.proof==NetworkClosureProof::ExactFullRankObstruction) {
            if (rank!=m || !cert.rigid_kernel.empty()) return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        const std::size_t nullity=m-rank;

        if (nullity==1) {
            if (cert.rigid_kernel.size()!=m) return invalid_certificate(cert);
            const auto kernel=decode_vector(sys.field,cert.rigid_kernel);
            if (!matrix_times_vector_zero(sys.field,sys.matrix,kernel))
                return invalid_certificate(cert);

            bool nonzero=false;
            for (const auto& x : kernel) if (!sys.field.is_zero(x)) nonzero=true;
            if (!nonzero) return invalid_certificate(cert);

            const auto signs=summarize_signs(sys.field,kernel);

            if (cert.status==NetworkClosureStatus::Closed &&
                cert.proof==NetworkClosureProof::CertifiedPositiveRigidKernel) {
                if (signs.indeterminate || signs.has_zero || signs.has_neg || !signs.has_pos)
                    return invalid_certificate(cert);
                return result_from_certificate(cert);
            }

            if (cert.status==NetworkClosureStatus::NotClosed &&
                cert.proof==NetworkClosureProof::CertifiedRigidSignObstruction) {
                if (signs.indeterminate) return invalid_certificate(cert);
                if (!(signs.has_zero || (signs.has_pos && signs.has_neg)))
                    return invalid_certificate(cert);
                return result_from_certificate(cert);
            }

            if (cert.status==NetworkClosureStatus::Indeterminate &&
                cert.proof==NetworkClosureProof::SignCertificationLimit) {
                if (!signs.indeterminate) return invalid_certificate(cert);
                return result_from_certificate(cert);
            }
            return invalid_certificate(cert);
        }

        if (nullity>1) {
            std::vector<Alg> recomputed_kernel;
            const auto cone=solve_general_positive_kernel(
                sys.field,rr,sys.matrix,m,&recomputed_kernel);

            if (cert.status==NetworkClosureStatus::Closed &&
                cert.proof==NetworkClosureProof::CertifiedGeneralPositiveKernel) {
                if (cone.state!=ConeState::Feasible ||
                    cert.rigid_kernel.size()!=m)
                    return invalid_certificate(cert);
                const auto kernel=decode_vector(sys.field,cert.rigid_kernel);
                if (!matrix_times_vector_zero(sys.field,sys.matrix,kernel))
                    return invalid_certificate(cert);
                const auto signs=summarize_signs(sys.field,kernel);
                if (signs.indeterminate || signs.has_zero || signs.has_neg || !signs.has_pos)
                    return invalid_certificate(cert);
                return result_from_certificate(cert);
            }

            if (cert.status==NetworkClosureStatus::NotClosed &&
                cert.proof==NetworkClosureProof::CertifiedGeneralConeObstruction) {
                if (cone.state!=ConeState::Infeasible || !cert.rigid_kernel.empty())
                    return invalid_certificate(cert);
                return result_from_certificate(cert);
            }

            if (cert.status==NetworkClosureStatus::Indeterminate &&
                cert.proof==NetworkClosureProof::EliminationComplexityLimit) {
                if (cone.state!=ConeState::ComplexityLimit || !cert.rigid_kernel.empty())
                    return invalid_certificate(cert);
                return result_from_certificate(cert);
            }

            if (cert.status==NetworkClosureStatus::Indeterminate &&
                cert.proof==NetworkClosureProof::SignCertificationLimit) {
                if (cone.state!=ConeState::Indeterminate)
                    return invalid_certificate(cert);
                return result_from_certificate(cert);
            }

            // Backward compatibility for older v1.0 certificates that stopped
            // at higher-dimensional kernel detection.
            if (cert.status==NetworkClosureStatus::Indeterminate &&
                cert.proof==NetworkClosureProof::HigherDimensionalKernel &&
                cert.rigid_kernel.empty())
                return result_from_certificate(cert);

            return invalid_certificate(cert);
        }

        return invalid_certificate(cert);
    } catch (...) {
        return invalid_certificate(cert);
    }
}

} // namespace pcg
