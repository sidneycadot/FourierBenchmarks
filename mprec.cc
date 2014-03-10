
#include <cstdio>
#include <iostream>
#include <mpfr.h>

using namespace std;

//MPFR_RNDN

int main()
{
    cout << "MPFR_PREC_MIN: " << MPFR_PREC_MIN << endl;
    cout << "MPFR_PREC_MAX: " << MPFR_PREC_MAX << endl;

    mpfr_t x;

    mpfr_init2(x, 1000);

    mpfr_const_pi(x, MPFR_RNDN);

    mpfr_printf("x: <%.1000Rf>\n", x);

    mpfr_clear(x);

    mpfr_free_cache();

    return 0;
}
