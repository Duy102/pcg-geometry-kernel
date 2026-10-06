#include <pcg/abcabc.hpp>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace pcg;

static PiRational parse_rational(const std::string& s) {
    auto slash=s.find('/');
    if (slash==std::string::npos) throw std::runtime_error("bad rational");
    return PiRational{BigInt{s.substr(0,slash)}, BigInt{s.substr(slash+1)}};
}
static Decision parse_decision(const std::string& s) {
    if(s=="REALIZABLE") return Decision::Realizable;
    if(s=="NOT_REALIZABLE") return Decision::NotRealizable;
    if(s=="INDETERMINATE") return Decision::Indeterminate;
    return Decision::Unsupported;
}
static int assurance_rank(ArithmeticAssurance a) {
    switch(a){case ArithmeticAssurance::Exact:return 2;case ArithmeticAssurance::CertifiedNumerical:return 1;case ArithmeticAssurance::None:return 0;}
    return -1;
}
static ArithmeticAssurance parse_assurance(const std::string& s) {
    if(s=="EXACT") return ArithmeticAssurance::Exact;
    if(s=="CERTIFIED_NUMERICAL") return ArithmeticAssurance::CertifiedNumerical;
    return ArithmeticAssurance::None;
}
int main(int argc,char**argv){
    if(argc!=2){std::cerr<<"corpus path required\n";return 2;}
    std::ifstream in(argv[1]);
    if(!in){std::cerr<<"cannot open corpus\n";return 2;}
    std::string line; int row=0,fail=0;
    while(std::getline(in,line)){
        if(line.empty()||line[0]=='#') continue;
        ++row;
        std::vector<std::string> f; std::stringstream ss(line); std::string x;
        while(std::getline(ss,x,',')) f.push_back(x);
        if(f.size()!=10){std::cerr<<"bad corpus row "<<row<<"\n";++fail;continue;}
        ABCABCInput input{{Turn{parse_rational(f[1])},Turn{parse_rational(f[2])},Turn{parse_rational(f[3])},Turn{parse_rational(f[4])},Turn{parse_rational(f[5])},Turn{parse_rational(f[6])}}};
        auto result=solve_abcabc(input);
        auto expected=parse_decision(f[7]);
        auto minimum=parse_assurance(f[8]);
        if(result.decision!=expected){std::cerr<<f[0]<<": decision mismatch\n";++fail;}
        if(expected!=Decision::Indeterminate && assurance_rank(result.assurance)<assurance_rank(minimum)){
            std::cerr<<f[0]<<": assurance below minimum\n";++fail;
        }
        auto verified=verify_abcabc_certificate(result.certificate);
        if(verified.decision!=result.decision){std::cerr<<f[0]<<": verifier mismatch\n";++fail;}
    }
    if(fail) return 1;
    std::cout<<"Golden corpus passed: "<<row<<" cases\n";
    return 0;
}
