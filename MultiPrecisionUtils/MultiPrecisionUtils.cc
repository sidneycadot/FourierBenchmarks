
////////////////////////////
// MultiPrecisionUtils.cc //
////////////////////////////

#include <cassert>

#include "MultiPrecisionUtils.h"

std::string to_string(const mpfr_t & x)
{
    char * str;
    int mpfr_asprintf_result = mpfr_asprintf(&str, "%.6Re", x);
    assert(mpfr_asprintf_result > 0);

    std::string r = str;

    mpfr_free_str(str);

    return r;
}

void gmp_randseed_string(gmp_randstate_t & state, const std::string & s)
{
    mpz_t seed;

    mpz_init(seed);

    for (unsigned i = 0; s[i] != '\0'; ++i)
    {
        const unsigned c = static_cast<unsigned char>(s[i]);

        mpz_mul_ui(seed, seed, 256);
        mpz_add_ui(seed, seed, c);
    }

    gmp_randseed(state, seed);

    mpz_clear(seed);
}
