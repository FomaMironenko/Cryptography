/* 
 * Decode pkcs v1.5 encoded message
 * clang++ -std=c++17 pkcs.cpp -lgmp -lgmpxx -o pkcs.o
 */

#include <iostream>
#include <algorithm>
#include <vector>

#include <gmp.h>
#include <gmpxx.h>


#define BYTE_SIZE 256

// RSA modulo and it's factorisation
const mpz_class N {
    "17976931348623159077293051907890247336179769789423065727343008115"
    "77326758055056206869853794492129829595855013875371640157101398586"
    "47833778606925583497541085196591615128057575940752635007475935288"
    "71082364994994077189561705436114947486504671101510156394068052754"
    "0071584560878577663743040086340742855278549092581"
};
const mpz_class p {
    "134078079299425970995740249982058461274793658205923933777235614437"
    "217640300736627688911116143623269986750405460943393208384195233759"
    "86027530441562135724301"
};
const mpz_class q = N / p;

// ciphertext
const mpz_class ct {
    "2209645186741038177630656113488341801741006978789283107173183914367"
    "6135600120538004282329650473509424343946219751512256465839967942889"
    "4607645420405815647489880137348641204523252293201764879166664029975"
    "0918872997169052608322206777160001932926087000957999372407745896777"
    "3697817571267229951148662959627934791540"
};
const mpz_class e = 65537;




// reduce pt by cropping the minor byte and return this byte as ASCII charcode
uint8_t bite_last_byte (
    mpz_class & pt
) {
    mpz_class r = pt % BYTE_SIZE;
    pt /= BYTE_SIZE;
    assert( r.fits_uint_p() );
    unsigned long r_ul = r.get_ui();
    assert( r_ul < BYTE_SIZE );
    return static_cast<uint8_t> ( r_ul );
}

std::vector<uint8_t> decode_PKCS (
    const mpz_class & ct,
    const mpz_class & p,
    const mpz_class & q,
    const mpz_class & pub_key
) {
    mpz_class N = p * q;
    mpz_class euler_f = (p - 1) * (q - 1); // euler function

    // compute the RSA secret key
    mpz_class sec_key{0};
    int ok = mpz_invert(
        sec_key.get_mpz_t(),
        pub_key.get_mpz_t(),
        euler_f.get_mpz_t()
    );
    assert(ok != 0 && "provided public key is not invertible");

    // compute the RSA plaintext
    mpz_class pt{0}; // plaintext
    mpz_powm(
        pt.get_mpz_t(),
        ct.get_mpz_t(),
        sec_key.get_mpz_t(),
        N.get_mpz_t()
    );

    // decode the plaintext
    // possible separators between the message and header
    const uint8_t separator1 = 0;
    const uint8_t separator2 = 255;
    std::vector<uint8_t> answer;
    uint8_t cur_byte = 0;
    while( pt != 0 ) {
        cur_byte = bite_last_byte(pt);
        if (cur_byte == separator1 || cur_byte == separator2) break;
        answer.push_back(cur_byte);
    }
    if (pt == 0) {
        std::cout << "WARNING: end reached without separator" << std::endl;
        return std::vector<uint8_t> (0);
    }

    // check PKCS consistency
    while( pt != 0 ) cur_byte = bite_last_byte(pt);
    if (cur_byte != 2 /*check last byte*/) {
        std::cout << "WARNING: provided ciphertext is not PKCS mod 2" << std::endl;
        return std::vector<uint8_t> (0);
    }

    std::reverse(answer.begin(), answer.end());
    return answer;
}


int main() {
    assert(N == p * q && "p does not divide N");

    std::vector<uint8_t> answer
        = decode_PKCS(ct, p, q, e);

    if (!answer.empty()) {
        std::cout 
            << "Ciphertext successfully decrypted. All checks passed. The plaintext is:" << std::endl 
            << answer.data() << std::endl;
    }

    return EXIT_SUCCESS;
}
