#include <pcg/abcabc.hpp>
#include <array>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace pcg;

namespace {
int failures = 0;

void check(bool cond, const std::string& name) {
    if (!cond) {
        std::cerr << "FAIL: " << name << "\n";
        ++failures;
    }
}

Decision parse_decision(const std::string& s) {
    if (s == "REALIZABLE") return Decision::Realizable;
    if (s == "NOT_REALIZABLE") return Decision::NotRealizable;
    throw std::runtime_error("bad decision in sign-rotation corpus");
}

ABCABCSignPattern parse_pattern(const std::string& s) {
    if (s.size() != 6) throw std::runtime_error("bad sign pattern length");
    ABCABCSignPattern out{};
    for (std::size_t i=0;i<6;++i) {
        if (s[i] == '+') out[i] = +1;
        else if (s[i] == '-') out[i] = -1;
        else throw std::runtime_error("bad sign character");
    }
    return out;
}

ABCABCSignPattern signs_of(const std::array<PiRational,6>& q) {
    ABCABCSignPattern out{};
    for (std::size_t i=0;i<6;++i)
        out[i] = q[i].value > Rational(BigInt{0}) ? +1 : -1;
    return out;
}

std::array<PiRational,6> rotate_q(const std::array<PiRational,6>& q, std::size_t shift) {
    std::array<PiRational,6> out{};
    for (std::size_t i=0;i<6;++i) out[i] = q[(i+shift)%6];
    return out;
}

std::array<PiRational,6> negate_q(const std::array<PiRational,6>& q) {
    std::array<PiRational,6> out{};
    for (std::size_t i=0;i<6;++i) out[i] = -q[i];
    return out;
}

ABCABCInput make_input(const std::array<PiRational,6>& q) {
    return ABCABCInput{{Turn{q[0]},Turn{q[1]},Turn{q[2]},Turn{q[3]},Turn{q[4]},Turn{q[5]}}};
}

void check_rotation(const std::array<PiRational,6>& q, int rotation, const std::string& label) {
    Rational sum(BigInt{0});
    for (const auto& x : q) sum += x.value;
    check(sum == Rational(BigInt{2*rotation}), label + " rotation sum");
}

void add_witness_orbit(const std::array<PiRational,6>& representative,
                       int rotation,
                       std::set<std::string>& seen,
                       int& witness_count,
                       const std::string& label) {
    for (std::size_t shift=0;shift<6;++shift) {
        auto q = rotate_q(representative,shift);
        auto input = make_input(q);
        const auto key = canonicalize_abcabc_input(input);
        if (!seen.insert(key).second) continue;

        ++witness_count;
        check_rotation(q,rotation,label);
        check(classify_abcabc_sign_rotation(signs_of(q),rotation) == Decision::Realizable,
              label + " classifier accepts theorem witness");
        const auto solved = solve_abcabc(input);
        check(solved.decision == Decision::Realizable,
              label + " concrete witness solves REALIZABLE");
        check(verify_abcabc_certificate(solved.certificate).decision == Decision::Realizable,
              label + " concrete witness certificate verifies");
    }
}
}

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "sign-rotation corpus path required\n";
        return 2;
    }

    // 1) Exhaustive theorem-level corpus: all 64 signs x rotations {-2,0,2}.
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "cannot open sign-rotation corpus\n";
        return 2;
    }

    int rows=0;
    int total_realizable=0;
    int realizable_minus2=0, realizable_zero=0, realizable_plus2=0;
    std::string line;
    while (std::getline(in,line)) {
        if (line.empty() || line[0]=='#') continue;
        std::stringstream ss(line);
        std::string pattern_s, rotation_s, expected_s;
        if (!std::getline(ss,pattern_s,',') ||
            !std::getline(ss,rotation_s,',') ||
            !std::getline(ss,expected_s,',')) {
            std::cerr << "bad corpus row\n";
            return 2;
        }
        ++rows;
        const int rotation = std::stoi(rotation_s);
        const Decision expected = parse_decision(expected_s);
        const Decision actual = classify_abcabc_sign_rotation(parse_pattern(pattern_s),rotation);
        check(actual == expected, "192-class corpus row " + pattern_s + "," + rotation_s);
        if (actual == Decision::Realizable) {
            ++total_realizable;
            if (rotation == -2) ++realizable_minus2;
            if (rotation == 0) ++realizable_zero;
            if (rotation == 2) ++realizable_plus2;
        }
    }
    check(rows == 192, "sign-rotation corpus has exactly 192 classes");
    check(total_realizable == 36, "classification has exactly 36 realizable classes");
    check(realizable_minus2 == 15, "r=-2 has exactly 15 realizable classes");
    check(realizable_zero == 6, "r=0 has exactly 6 realizable classes");
    check(realizable_plus2 == 15, "r=2 has exactly 15 realizable classes");

    // Other integer rotations are rejected by the theorem's Seifert filter.
    check(classify_abcabc_sign_rotation(parse_pattern("++++++"),1) == Decision::NotRealizable,
          "unsupported theorem rotation r=1 is impossible");
    check(classify_abcabc_sign_rotation(parse_pattern("++++++"),3) == Decision::NotRealizable,
          "unsupported theorem rotation r=3 is impossible");

    bool invalid_sign_rejected=false;
    try {
        (void)classify_abcabc_sign_rotation(ABCABCSignPattern{{1,1,1,1,1,0}},2);
    } catch (const std::invalid_argument&) {
        invalid_sign_rejected=true;
    }
    check(invalid_sign_rejected, "invalid sign alphabet rejected");

    // 2) Concrete fixed-turn witnesses from v0.4.
    // Manuscript uses x=tau/(2*pi); q below is tau/pi = 2x.
    const std::array<PiRational,6> r0{{
        PiRational{1,2}, PiRational{7,4}, PiRational{1,2},
        PiRational{-1,1}, PiRational{-3,4}, PiRational{-1,1}
    }};
    const std::array<std::array<PiRational,6>,4> r2{{
        std::array<PiRational,6>{{PiRational{1,3},PiRational{1,1},PiRational{1,3},PiRational{1,1},PiRational{1,3},PiRational{1,1}}},
        std::array<PiRational,6>{{PiRational{5,4},PiRational{1,2},PiRational{3,4},PiRational{1,4},PiRational{3,2},PiRational{-1,4}}},
        std::array<PiRational,6>{{PiRational{5,4},PiRational{1,4},PiRational{5,4},PiRational{-1,4},PiRational{7,4},PiRational{-1,4}}},
        std::array<PiRational,6>{{PiRational{3,2},PiRational{-1,6},PiRational{3,2},PiRational{-1,6},PiRational{3,2},PiRational{-1,6}}}
    }};

    std::set<std::string> seen_plus2;
    int plus2_witnesses=0;
    for (std::size_t i=0;i<r2.size();++i)
        add_witness_orbit(r2[i],2,seen_plus2,plus2_witnesses,"r=2 orbit " + std::to_string(i+1));
    check(plus2_witnesses == 15, "r=2 witnesses cover exactly 15 sign classes");

    std::set<std::string> seen_minus2;
    int minus2_witnesses=0;
    for (std::size_t i=0;i<r2.size();++i)
        add_witness_orbit(negate_q(r2[i]),-2,seen_minus2,minus2_witnesses,"r=-2 reflected orbit " + std::to_string(i+1));
    check(minus2_witnesses == 15, "r=-2 witnesses cover exactly 15 sign classes");

    std::set<std::string> seen_zero;
    int zero_witnesses=0;
    add_witness_orbit(r0,0,seen_zero,zero_witnesses,"r=0 orbit");
    check(zero_witnesses == 6, "r=0 witnesses cover exactly 6 sign classes");

    check(plus2_witnesses + minus2_witnesses + zero_witnesses == 36,
          "concrete theorem witnesses cover all 36 realizable sign-rotation classes");

    if (failures) {
        std::cerr << failures << " sign-rotation test(s) failed\n";
        return 1;
    }
    std::cout << "ABCABC sign-rotation corpus passed: 192 classes; 36 concrete realizable witnesses verified\n";
    return 0;
}
