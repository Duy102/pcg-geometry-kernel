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
const char* kSchemaVersion = "1.1";
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

enum class FeasibilityStatus { Feasible, Infeasible, Indeterminate };

struct FeasibilityResult {
    FeasibilityStatus status{FeasibilityStatus::Indeterminate};
    std::vector<Alg> variables;
};

FeasibilityResult nonnegative_equality_feasibility(
        const CyclotomicField& field,
        const std::vector<std::vector<Alg>>& input_a,
        const std::vector<Alg>& input_b,
        std::size_t variable_count) {
    if (input_a.size()!=input_b.size())
        throw std::invalid_argument("equality system row/rhs mismatch");
    for (const auto& row : input_a)
        if (row.size()!=variable_count)
            throw std::invalid_argument("equality system column mismatch");

    std::vector<std::vector<Alg>> a;
    std::vector<Alg> rhs;
    a.reserve(input_a.size());
    rhs.reserve(input_b.size());

    for (std::size_t i=0;i<input_a.size();++i) {
        bool row_zero=true;
        for (const auto& x : input_a[i])
            if (!field.is_zero(x)) { row_zero=false; break; }

        if (row_zero) {
            if (!field.is_zero(input_b[i]))
                return {FeasibilityStatus::Infeasible,{}};
            continue;
        }

        auto row=input_a[i];
        Alg b=input_b[i];
        const auto sign=certified_real_sign(field,b);
        if (sign==SignCert::Indeterminate)
            return {FeasibilityStatus::Indeterminate,{}};
        if (sign==SignCert::Negative) {
            for (auto& x : row) x=field.neg(x);
            b=field.neg(b);
        }
        a.push_back(std::move(row));
        rhs.push_back(std::move(b));
    }

    const std::size_t rows=a.size();
    if (rows==0)
        return {FeasibilityStatus::Feasible,
                std::vector<Alg>(variable_count,field.zero())};

    const std::size_t total_vars=variable_count+rows;
    std::vector<std::vector<Alg>> tab(
        rows,std::vector<Alg>(total_vars,field.zero()));
    std::vector<std::size_t> basis(rows);

    for (std::size_t i=0;i<rows;++i) {
        for (std::size_t j=0;j<variable_count;++j)
            tab[i][j]=a[i][j];
        tab[i][variable_count+i]=field.one();
        basis[i]=variable_count+i;
    }

    auto cost_is_one=[&](std::size_t var) {
        return var>=variable_count;
    };

    const std::size_t iteration_limit=100000;
    for (std::size_t iteration=0;iteration<iteration_limit;++iteration) {
        std::vector<bool> is_basic(total_vars,false);
        for (auto b : basis) is_basic[b]=true;

        std::size_t enter=total_vars;
        for (std::size_t j=0;j<total_vars;++j) {
            if (is_basic[j]) continue;

            Alg reduced=cost_is_one(j)?field.one():field.zero();
            for (std::size_t i=0;i<rows;++i) {
                if (cost_is_one(basis[i]))
                    reduced=field.sub(reduced,tab[i][j]);
            }

            const auto sign=certified_real_sign(field,reduced);
            if (sign==SignCert::Indeterminate)
                return {FeasibilityStatus::Indeterminate,{}};
            if (sign==SignCert::Negative) {
                enter=j; // Bland: first eligible index.
                break;
            }
        }

        if (enter==total_vars) {
            Alg objective=field.zero();
            for (std::size_t i=0;i<rows;++i)
                if (cost_is_one(basis[i]))
                    objective=field.add(objective,rhs[i]);

            if (field.is_zero(objective)) {
                std::vector<Alg> x(variable_count,field.zero());
                for (std::size_t i=0;i<rows;++i)
                    if (basis[i]<variable_count)
                        x[basis[i]]=rhs[i];

                for (const auto& v : x) {
                    const auto sign=certified_real_sign(field,v);
                    if (sign==SignCert::Negative)
                        return {FeasibilityStatus::Indeterminate,{}};
                    if (sign==SignCert::Indeterminate)
                        return {FeasibilityStatus::Indeterminate,{}};
                }
                return {FeasibilityStatus::Feasible,std::move(x)};
            }

            const auto sign=certified_real_sign(field,objective);
            if (sign==SignCert::Positive)
                return {FeasibilityStatus::Infeasible,{}};
            return {FeasibilityStatus::Indeterminate,{}};
        }

        std::size_t leave=rows;
        Alg best_ratio=field.zero();
        bool have_ratio=false;
        for (std::size_t i=0;i<rows;++i) {
            const auto coeff_sign=certified_real_sign(field,tab[i][enter]);
            if (coeff_sign==SignCert::Indeterminate)
                return {FeasibilityStatus::Indeterminate,{}};
            if (coeff_sign!=SignCert::Positive) continue;

            const Alg ratio=field.divide(rhs[i],tab[i][enter]);
            if (!have_ratio) {
                leave=i;
                best_ratio=ratio;
                have_ratio=true;
                continue;
            }

            const auto cmp=certified_real_sign(field,field.sub(ratio,best_ratio));
            if (cmp==SignCert::Indeterminate)
                return {FeasibilityStatus::Indeterminate,{}};
            if (cmp==SignCert::Negative ||
                (cmp==SignCert::Zero && basis[i]<basis[leave])) {
                leave=i;
                best_ratio=ratio;
            }
        }

        if (!have_ratio)
            return {FeasibilityStatus::Indeterminate,{}};

        const Alg pivot=tab[leave][enter];
        const Alg inv=field.inverse(pivot);
        for (std::size_t j=0;j<total_vars;++j)
            tab[leave][j]=field.mul(tab[leave][j],inv);
        rhs[leave]=field.mul(rhs[leave],inv);

        for (std::size_t i=0;i<rows;++i) {
            if (i==leave || field.is_zero(tab[i][enter])) continue;
            const Alg factor=tab[i][enter];
            for (std::size_t j=0;j<total_vars;++j)
                tab[i][j]=field.sub(tab[i][j],field.mul(factor,tab[leave][j]));
            rhs[i]=field.sub(rhs[i],field.mul(factor,rhs[leave]));
        }
        basis[leave]=enter;
    }

    return {FeasibilityStatus::Indeterminate,{}};
}

FeasibilityResult find_positive_kernel(const CyclotomicField& field,
                                       const std::vector<std::vector<Alg>>& matrix,
                                       std::size_t columns) {
    std::vector<Alg> rhs(matrix.size(),field.zero());
    for (std::size_t i=0;i<matrix.size();++i) {
        Alg row_sum=field.zero();
        for (std::size_t j=0;j<columns;++j)
            row_sum=field.add(row_sum,matrix[i][j]);
        rhs[i]=field.neg(row_sum);
    }

    auto phase=nonnegative_equality_feasibility(field,matrix,rhs,columns);
    if (phase.status!=FeasibilityStatus::Feasible)
        return phase;

    for (auto& x : phase.variables)
        x=field.add(x,field.one()); // c = 1 + z, hence c is strictly positive.
    return phase;
}

FeasibilityResult find_stiemke_dual(const CyclotomicField& field,
                                    const std::vector<std::vector<Alg>>& matrix,
                                    std::size_t columns) {
    const std::size_t rows=matrix.size();
    const std::size_t vars=2*rows+columns;
    std::vector<std::vector<Alg>> a(
        columns+1,std::vector<Alg>(vars,field.zero()));
    std::vector<Alg> b(columns+1,field.zero());

    // z = N^T(y+ - y-) >= 0.
    for (std::size_t e=0;e<columns;++e) {
        for (std::size_t i=0;i<rows;++i) {
            a[e][i]=field.neg(matrix[i][e]);
            a[e][rows+i]=matrix[i][e];
        }
        a[e][2*rows+e]=field.one();
    }

    // Normalize the nonzero nonnegative stress: sum(z)=1.
    for (std::size_t e=0;e<columns;++e)
        a[columns][2*rows+e]=field.one();
    b[columns]=field.one();

    auto phase=nonnegative_equality_feasibility(field,a,b,vars);
    if (phase.status!=FeasibilityStatus::Feasible)
        return phase;

    std::vector<Alg> y(rows,field.zero());
    for (std::size_t i=0;i<rows;++i)
        y[i]=field.sub(phase.variables[i],phase.variables[rows+i]);
    return {FeasibilityStatus::Feasible,std::move(y)};
}

bool positive_kernel_valid(const CyclotomicField& field,
                           const std::vector<std::vector<Alg>>& matrix,
                           const std::vector<Alg>& c,
                           bool* comparison_indeterminate=nullptr) {
    if (comparison_indeterminate) *comparison_indeterminate=false;
    if (matrix.empty()) {
        for (const auto& x : c) {
            const auto sign=certified_real_sign(field,x);
            if (sign==SignCert::Indeterminate) {
                if (comparison_indeterminate) *comparison_indeterminate=true;
                return false;
            }
            if (sign!=SignCert::Positive) return false;
        }
        return true;
    }
    if (matrix[0].size()!=c.size()) return false;
    if (!matrix_times_vector_zero(field,matrix,c)) return false;

    for (const auto& x : c) {
        const auto sign=certified_real_sign(field,x);
        if (sign==SignCert::Indeterminate) {
            if (comparison_indeterminate) *comparison_indeterminate=true;
            return false;
        }
        if (sign!=SignCert::Positive) return false;
    }
    return true;
}

bool stiemke_dual_valid(const CyclotomicField& field,
                        const std::vector<std::vector<Alg>>& matrix,
                        const std::vector<Alg>& y,
                        bool* comparison_indeterminate=nullptr) {
    if (comparison_indeterminate) *comparison_indeterminate=false;
    const std::size_t rows=matrix.size();
    const std::size_t cols=rows?matrix[0].size():0;
    if (y.size()!=rows || cols==0) return false;

    for (const auto& v : y)
        if (!field.is_real(v)) return false;

    bool positive=false;
    for (std::size_t e=0;e<cols;++e) {
        Alg z=field.zero();
        for (std::size_t i=0;i<rows;++i)
            z=field.add(z,field.mul(matrix[i][e],y[i]));

        const auto sign=certified_real_sign(field,z);
        if (sign==SignCert::Indeterminate) {
            if (comparison_indeterminate) *comparison_indeterminate=true;
            return false;
        }
        if (sign==SignCert::Negative) return false;
        if (sign==SignCert::Positive) positive=true;
    }
    return positive;
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
    else if (cert.proof==NetworkClosureProof::CertifiedPositiveKernel ||
             cert.proof==NetworkClosureProof::CertifiedStiemkeObstruction)
        out.assurance=ArithmeticAssurance::CertifiedNumerical;
    else
        out.assurance=ArithmeticAssurance::None;
    out.termination=cert.status==NetworkClosureStatus::Indeterminate
        ? (cert.proof==NetworkClosureProof::CyclotomicOrderLimit ||
           cert.proof==NetworkClosureProof::SignCertificationLimit
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

std::pair<Alg,Alg> exact_unit_direction(const CyclotomicField& field,
                                        PiRational phase) {
    if (field.order%4U!=0U)
        throw std::runtime_error("cyclotomic order is incompatible with pi-phase encoding");
    const BigInt L=BigInt{field.order/4U};
    const BigInt scaled_num=phase.value.numerator()*L;
    const BigInt den=phase.value.denominator();
    if (scaled_num%den!=0)
        throw std::runtime_error("phase is not representable in the active cyclotomic field");

    BigInt exp_big=BigInt{2}*(scaled_num/den);
    exp_big%=field.order;
    if (exp_big<0) exp_big+=field.order;
    const unsigned exp=exp_big.convert_to<unsigned>();
    const unsigned negexp=(field.order-exp)%field.order;

    const Alg z=field.monomial(exp);
    const Alg zbar=field.monomial(negexp);
    const Rational half(BigInt{1},BigInt{2});
    const Alg imag_unit=field.monomial(field.order/4U);
    const Alg cosine=field.scale(field.add(z,zbar),half);
    const Alg sine=field.scale(
        field.mul(field.sub(z,zbar),field.neg(imag_unit)),
        half);

    if (!field.is_real(cosine) || !field.is_real(sine))
        throw std::runtime_error("exact unit direction is not real");
    return {cosine,sine};
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
    auto emit_vector=[&](std::ostringstream& os,
                         const std::vector<AlgebraicChordValue>& values) {
        os << "[";
        for (std::size_t i=0;i<values.size();++i) {
            if (i) os << ",";
            os << "[";
            for (std::size_t j=0;j<values[i].coefficients.size();++j) {
                if (j) os << ",";
                os << "\"" << rational_text(values[i].coefficients[j]) << "\"";
            }
            os << "]";
        }
        os << "]";
    };

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
       << ",\"positive_kernel\":";
    emit_vector(os,c.positive_kernel);
    os << ",\"stiemke_dual\":";
    emit_vector(os,c.stiemke_dual);
    os << "}";
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
        } else {
            auto positive=find_positive_kernel(sys.field,sys.matrix,m);
            if (positive.status==FeasibilityStatus::Indeterminate) {
                cert.status=NetworkClosureStatus::Indeterminate;
                cert.proof=NetworkClosureProof::SignCertificationLimit;
                return result_from_certificate(cert);
            }

            if (positive.status==FeasibilityStatus::Feasible) {
                bool sign_limit=false;
                if (!positive_kernel_valid(sys.field,sys.matrix,positive.variables,&sign_limit)) {
                    cert.status=NetworkClosureStatus::Indeterminate;
                    cert.proof=sign_limit
                        ? NetworkClosureProof::SignCertificationLimit
                        : NetworkClosureProof::BackendFailure;
                    return result_from_certificate(cert);
                }
                cert.status=NetworkClosureStatus::Closed;
                cert.proof=NetworkClosureProof::CertifiedPositiveKernel;
                cert.positive_kernel=encode_vector(positive.variables);
            } else {
                auto dual=find_stiemke_dual(sys.field,sys.matrix,m);
                if (dual.status!=FeasibilityStatus::Feasible) {
                    cert.status=NetworkClosureStatus::Indeterminate;
                    cert.proof=dual.status==FeasibilityStatus::Indeterminate
                        ? NetworkClosureProof::SignCertificationLimit
                        : NetworkClosureProof::BackendFailure;
                    return result_from_certificate(cert);
                }

                bool sign_limit=false;
                if (!stiemke_dual_valid(sys.field,sys.matrix,dual.variables,&sign_limit)) {
                    cert.status=NetworkClosureStatus::Indeterminate;
                    cert.proof=sign_limit
                        ? NetworkClosureProof::SignCertificationLimit
                        : NetworkClosureProof::BackendFailure;
                    return result_from_certificate(cert);
                }
                cert.status=NetworkClosureStatus::NotClosed;
                cert.proof=NetworkClosureProof::CertifiedStiemkeObstruction;
                cert.stiemke_dual=encode_vector(dual.variables);
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
                cert.positive_kernel.empty() &&
                cert.stiemke_dual.empty())
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
            if (rank!=m ||
                !cert.positive_kernel.empty() ||
                !cert.stiemke_dual.empty())
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        if (cert.status==NetworkClosureStatus::Closed &&
            cert.proof==NetworkClosureProof::CertifiedPositiveKernel) {
            if (rank>=m || cert.positive_kernel.size()!=m ||
                !cert.stiemke_dual.empty())
                return invalid_certificate(cert);
            const auto c=decode_vector(sys.field,cert.positive_kernel);
            bool sign_limit=false;
            if (!positive_kernel_valid(sys.field,sys.matrix,c,&sign_limit))
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        if (cert.status==NetworkClosureStatus::NotClosed &&
            cert.proof==NetworkClosureProof::CertifiedStiemkeObstruction) {
            if (rank>=m || !cert.positive_kernel.empty() ||
                cert.stiemke_dual.size()!=sys.matrix.size())
                return invalid_certificate(cert);
            const auto y=decode_vector(sys.field,cert.stiemke_dual);
            bool sign_limit=false;
            if (!stiemke_dual_valid(sys.field,sys.matrix,y,&sign_limit))
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        if (cert.status==NetworkClosureStatus::Indeterminate &&
            cert.proof==NetworkClosureProof::SignCertificationLimit) {
            if (rank>=m ||
                !cert.positive_kernel.empty() ||
                !cert.stiemke_dual.empty())
                return invalid_certificate(cert);
            return result_from_certificate(cert);
        }

        return invalid_certificate(cert);
    } catch (...) {
        return invalid_certificate(cert);
    }
}


CertifiedScalar certified_algebraic_real(
    unsigned cyclotomic_order,
    const AlgebraicChordValue& value) {
    if (cyclotomic_order==0)
        throw std::invalid_argument("cyclotomic order must be positive");
    CyclotomicField field(cyclotomic_order);
    if (value.coefficients.size()!=field.degree)
        throw std::invalid_argument("algebraic coefficient count does not match cyclotomic degree");

    Alg a{value.coefficients};
    if (!field.is_real(a))
        throw std::invalid_argument("algebraic value is not exactly real");

    detail::Interval acc(0.0);
    for (std::size_t j=0;j<field.degree;++j) {
        if (a.c[j]==Rational(BigInt{0})) continue;
        const PiRational phase{BigInt{2}*BigInt{j},BigInt{cyclotomic_order}};
        acc += rational_interval_local(a.c[j]) * certified_cos_pi(phase).interval();
    }
    return CertifiedScalar(acc);
}

ProjectivelyRigidEmbeddingResult reconstruct_positive_kernel_embedding(
    const NetworkClosureCertificate& cert) {
    ProjectivelyRigidEmbeddingResult out;

    try {
        const auto verified=verify_network_closure_certificate(cert);
        if (verified.termination==TerminationReason::BackendFailure) {
            out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
            return out;
        }
        if (verified.status==NetworkClosureStatus::Indeterminate) {
            out.status=ProjectivelyRigidEmbeddingStatus::Indeterminate;
            return out;
        }
        if (verified.status!=NetworkClosureStatus::Closed) {
            out.status=ProjectivelyRigidEmbeddingStatus::NotClosed;
            return out;
        }

        const std::size_t m=cert.input.turns.size();

        auto sys=build_exact_system(cert.input);
        if (!sys.supported || sys.order!=cert.cyclotomic_order) {
            out.status=ProjectivelyRigidEmbeddingStatus::Indeterminate;
            return out;
        }

        auto c=decode_vector(sys.field,cert.positive_kernel);
        Alg total=sys.field.zero();
        for (const auto& x : c) total=sys.field.add(total,x);
        if (sys.field.is_zero(total)) {
            out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
            return out;
        }
        for (auto& x : c) x=sys.field.divide(x,total);

        std::vector<Alg> px(cert.input.vertex_count,sys.field.zero());
        std::vector<Alg> py(cert.input.vertex_count,sys.field.zero());
        std::vector<bool> assigned(cert.input.vertex_count,false);

        Alg curx=sys.field.zero();
        Alg cury=sys.field.zero();
        const std::size_t root=cert.input.trace_vertices.front();
        assigned[root]=true;

        for (std::size_t e=0;e<m;++e) {
            const std::size_t tail=cert.input.trace_vertices[e];
            const std::size_t head=cert.input.trace_vertices[(e+1)%m];

            if (!assigned[tail]) {
                px[tail]=curx;
                py[tail]=cury;
                assigned[tail]=true;
            } else if (!sys.field.equal(px[tail],curx) ||
                       !sys.field.equal(py[tail],cury)) {
                out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
                return out;
            }

            const auto [ux,uy]=exact_unit_direction(
                sys.field,sys.topology.chord_phase_pi[e]);
            const Alg nextx=sys.field.add(curx,sys.field.mul(c[e],ux));
            const Alg nexty=sys.field.add(cury,sys.field.mul(c[e],uy));

            if (!assigned[head]) {
                px[head]=nextx;
                py[head]=nexty;
                assigned[head]=true;
            } else if (!sys.field.equal(px[head],nextx) ||
                       !sys.field.equal(py[head],nexty)) {
                out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
                return out;
            }

            curx=nextx;
            cury=nexty;
        }

        if (!sys.field.equal(curx,px[root]) ||
            !sys.field.equal(cury,py[root])) {
            out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
            return out;
        }

        for (bool x : assigned) {
            if (!x) {
                out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
                return out;
            }
        }

        out.embedding.cyclotomic_order=sys.order;
        out.embedding.vertices.reserve(cert.input.vertex_count);
        for (std::size_t v=0;v<cert.input.vertex_count;++v)
            out.embedding.vertices.push_back(
                AlgebraicPoint2{AlgebraicChordValue{px[v].c},
                                AlgebraicChordValue{py[v].c}});

        out.embedding.chord_magnitudes=encode_vector(c);
        out.embedding.tangent_phase_pi.reserve(m);
        PiRational phase{0,1};
        for (std::size_t e=0;e<m;++e) {
            out.embedding.tangent_phase_pi.push_back(phase);
            phase=phase+cert.input.turns[e].pi;
        }

        out.status=ProjectivelyRigidEmbeddingStatus::Ready;
        return out;
    } catch (...) {
        out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
        return out;
    }
}


ProjectivelyRigidEmbeddingResult reconstruct_projectively_rigid_embedding(
    const NetworkClosureCertificate& cert) {
    try {
        const auto verified=verify_network_closure_certificate(cert);
        if (verified.termination==TerminationReason::BackendFailure) {
            ProjectivelyRigidEmbeddingResult out;
            out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
            return out;
        }
        if (verified.status==NetworkClosureStatus::Indeterminate) {
            ProjectivelyRigidEmbeddingResult out;
            out.status=ProjectivelyRigidEmbeddingStatus::Indeterminate;
            return out;
        }
        if (verified.status!=NetworkClosureStatus::Closed) {
            ProjectivelyRigidEmbeddingResult out;
            out.status=ProjectivelyRigidEmbeddingStatus::NotClosed;
            return out;
        }

        const std::size_t m=cert.input.turns.size();
        if (cert.exact_rank>=m || m-cert.exact_rank!=1) {
            ProjectivelyRigidEmbeddingResult out;
            out.status=ProjectivelyRigidEmbeddingStatus::NotProjectivelyRigid;
            return out;
        }
        return reconstruct_positive_kernel_embedding(cert);
    } catch (...) {
        ProjectivelyRigidEmbeddingResult out;
        out.status=ProjectivelyRigidEmbeddingStatus::InvalidCertificate;
        return out;
    }
}

} // namespace pcg
