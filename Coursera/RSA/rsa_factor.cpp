/* 
 * A toy program to compute RSA factorization.
 * Compile with
 * clang++ -std=c++17 rsa_factor.cpp -lgmp -lgmpxx -o rsa_factor.o
 */

#include <iostream>

#include <gmp.h>
#include <gmpxx.h>


const mpz_class N1 {
    "17976931348623159077293051907890247336179769789423065727343008115"
    "77326758055056206869853794492129829595855013875371640157101398586"
    "47833778606925583497541085196591615128057575940752635007475935288"
    "71082364994994077189561705436114947486504671101510156394068052754"
    "0071584560878577663743040086340742855278549092581"
};

const mpz_class N2 {
    "6484558428080716696628242653467722787263437207069762630604390703787"
    "9730861808111646271401527606141756919558732184025452065542490671989"
    "2428844841839353281972988531310511738648965962582821502504990264452"
    "1008852816733037111422964210278402893076574586452336833570778346897"
    "15838646088239640236866252211790085787877"
};

const mpz_class N3 {
    "72006226374735042527956443552558373833808445147399984182665305798191"
    "63556901883377904234086641876639384851752649940178970835240791356868"
    "77441155132015188279331812309091996246361896836573643119174094961348"
    "52463970788523879939683923036467667022162701835329944324119217381272"
    "9276147530748597302192751375739387929"
};


// There is no need to search more than 2^64 iterations. It's too long
using search_type = uint64_t;

constexpr
const search_type const_pow(search_type base, uint8_t a) 
noexcept {
    if (a == 0) {
        return 1;
    } else if (a % 2 == 0) {
        const search_type half = pow(base, a / 2);
        return half * half;
    } else {
        return base * pow(base, a-1);
    }
}

mpz_class smallest_divisor(mpz_class N, mpz_class divisor) {
    mpz_class divisor1 = N / divisor;
    return (divisor < divisor1) ? divisor : divisor1;
}


// Let N = p*q -- a RSA modulo
// We assume that |a*p - b*q| < c * N ^ 0.25
// Than there exists an elegant way of factoring N
mpz_class factor_with_hint(
    const mpz_class &N, 
    uint8_t a = 1, uint8_t b = 1,
    search_type max_search = 10
) {
    // check dummy case
    if (a != 1 && N % a == 0) return smallest_divisor(N, a);
    if (b != 1 && N % b == 0) return smallest_divisor(N, b);

    mpz_class N1 = N * a * b * 4;
    mpz_class A{0};
    mpz_sqrt(A.get_mpz_t(), N1.get_mpz_t());

    std::cout   << "Searching a divisor with " << max_search 
                << " iterations ..." << std::endl;
    mpz_class x{0};
    mpz_class d{0};
    for (search_type i = 0; i < max_search; i++, A += 1) {
        x = A*A - N1;
        if (x < 0) continue; // possible only at the first iteration
        mpz_sqrt(x.get_mpz_t(), x.get_mpz_t());
        d = A+x;
        if (N1 % d == 0) {
            if (d % a == 0) d /= a;
            if (d % b == 0) d /= b;
            if (d % 4 == 0) d /= 4;
            if (d % 2 == 0) d /= 2;
            std::cout << "DONE" << std::endl;
            return smallest_divisor(N, d);
        }
    }
    std::cout << "Failed to find a divisor" << std::endl;
    return 0;
}


void do_work(
    const mpz_class &N, 
    uint8_t a = 1, uint8_t b = 1,
    search_type max_search = 10
) {
    static unsigned counter = 1;
    std::cout << "==== TASK" << counter++ << " ====" << std::endl;
    mpz_class A = factor_with_hint(N, a, b, max_search);
    if (A == 0 || N % A != 0) {
        std::cout << "Didn't manage to factor" << std::endl;
    } else {
        std::cout << "Factofization computed. All check passed." << std::endl;
        std::cout << "The factor is: \n" << A << std::endl;
    }
    std::cout << std::endl;
}

int main() {

    do_work(N1);

    do_work(N2, 1, 1, const_pow(2, 20));

    do_work(N3, 3, 2);

    return EXIT_SUCCESS;
}
 