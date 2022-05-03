/*
 * A toy programm to compute discrete log with some assumptions
 * coplie with
 * clang++ -std=c++17 main.cpp -lgmp -lgmpxx -o main.o
 */

#include <iostream>
#include <map>

#include <gmp.h>
#include <gmpxx.h>


const mpz_class p(
    "1340780792994259709957402499820584612747"
    "9365820592393377723561443721764030073546"
    "9768018742981669034276900318581864860508"
    "53753882811946569946433649006084171"
);

const mpz_class g(
    "1171782988036620700951611759633536708855"
    "8084999998952205599979459063929499736583"
    "7466705721764714603129285948296754282794"
    "66566527115212748467589894601965568"
);

const mpz_class h(
    "3239475104050450443565264378728065788649"
    "0975209524495278347924529719819761432925"
    "5807385693795855318053287892800149470609"
    "7394108577585732452307673444020333"
);



using power_type = uint32_t;

constexpr
const power_type const_pow(power_type base, uint8_t a) 
noexcept {
    if (a == 0) {
        return 1;
    } else if (a % 2 == 0) {
        const power_type half = pow(base, a / 2);
        return half * half;
    } else {
        return base * pow(base, a-1);
    }
}


class CyclicGroup {
public:
    const mpz_class mod;
    const mpz_class g;
    static const power_type B;

    CyclicGroup(const mpz_class &mod, const mpz_class &generator):
        mod(mod), g(generator)
    { }

    mpz_class pow(const mpz_class &x) const {
        mpz_class result{0};
        mpz_powm(
            result.get_mpz_t(), // output
            g.get_mpz_t(),      // base
            x.get_mpz_t(),      // power
            mod.get_mpz_t()     // modulo
        );
        return result;
    }

    /* Compute discrete log in this cyclic group,
     * assuming that the log is within 0 and 2^40
     * 
     * Use the following idea:
     * if x < 2^40, then there exist a unique x0, x1 < 2^20 s.t.
     * x = x0*B + x1  where B = 2^20
     * 
     * We derive:
     * h / g^x1 == (g^B)^x0
     * and use 'meet in the middle' attack 
     */
    mpz_class adapted_log(const mpz_class &h) const {

        std::map<mpz_class, power_type> LUT;

        // compute lhs lookup table
        mpz_class lhs = h;
        mpz_class g_inv{0};
        int ok = mpz_invert(
            g_inv.get_mpz_t(),
            g.get_mpz_t(),
            mod.get_mpz_t()
        );
        assert(ok != 0 && "generator must be invertible");
        std::cout << "Computing lhs lookup table ..." << std::endl;
        for (power_type x1 = 0; x1 < B; x1++, lhs *= g_inv) {
            lhs %= mod;
            LUT[lhs] = x1;
        }
        std::cout << "DONE" << std::endl;

        // compute rhs
        const mpz_class base = pow(B);
        mpz_class rhs = 1;
        decltype(LUT.find(rhs)) cached_lhs;
        std::cout << "Computing rhs ..." << std::endl;
        for (power_type x0 = 0; x0 < B; x0++, rhs *= base) {
            rhs %= mod;
            cached_lhs = LUT.find(rhs);
            if (cached_lhs != LUT.end()) {
                power_type x1 = cached_lhs->second;
                std::cout << "DONE" << std::endl;
                return mpz_class{x0} * mpz_class{B} + mpz_class{x1};
            }
        }

        std::cout << "Didn't manage to find a discrete log for " << h << std::endl;
        return mpz_class{0};
    }
};

const power_type CyclicGroup::B = const_pow(2, 20);



int main() {

    const CyclicGroup group(p, g);

    mpz_class x = group.adapted_log(h);

    mpz_class h1 = group.pow(x);

    if (h1 != h) {
        std::cout << "Didn't manage to compute discrete log";
    } else {
        std::cout << "Discrete log computed. All checks passed" << std::endl;
        std::cout << "The valid answer is: " << x << std::endl;
    }
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
