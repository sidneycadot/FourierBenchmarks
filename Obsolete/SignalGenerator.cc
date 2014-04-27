
#include <iostream>
#include <cmath>
#include <vector>
#include <mpfr.h>
#include <mpc.h>
#include "print_collection.h"

using namespace std;

const mpfr_rnd_t default_rnd = MPFR_RNDN;
const mpfr_prec_t default_prec = 4000; // will suffice for 1000 decimal digits of precision.

vector<mpfr_t *> make_real_sinusoid(const mpfr_t & fs, const unsigned & n, const mpfr_t & f, const mpfr_prec_t & prec, const mpfr_rnd_t & rnd)
{
    vector<mpfr_t *> sigvec;
    sigvec.reserve(n);

    // calculate sin(2 * pi * f * i / fs);

    mpfr_t coeff;
    mpfr_init2(coeff, prec);

    mpfr_const_pi(coeff, rnd);
    mpfr_mul_ui(coeff, coeff, 2, rnd);
    mpfr_mul(coeff, coeff, f, rnd);
    mpfr_div(coeff, coeff, fs, rnd);

    //

    mpfr_t x;
    mpfr_init2(x, prec);

    for (unsigned i = 0; i < n; ++i)
    {
        mpfr_mul_ui(x, coeff, i, rnd);

        mpfr_t * sig = new mpfr_t[1];

        mpfr_init2(sig[0], prec);

        mpfr_sin(sig[0], x, rnd);

        sigvec.push_back(sig);
    }

    mpfr_clear(coeff);
    mpfr_clear(x);

    return sigvec;
}

vector<mpc_t *> make_complex_sinusoid(const mpfr_t & fs, const unsigned & n, const mpfr_t & f, const mpfr_prec_t & prec, const mpfr_rnd_t & rnd)
{
    vector<mpc_t *> sigvec;
    sigvec.reserve(n);

    // calculate sin(2 * pi * f * i / fs);

    mpfr_t coeff;
    mpfr_init2(coeff, prec);

    mpfr_const_pi(coeff, rnd);
    mpfr_mul_ui(coeff, coeff, 2, rnd);
    mpfr_mul(coeff, coeff, f, rnd);
    mpfr_div(coeff, coeff, fs, rnd);

    //

    mpfr_t x;
    mpfr_init2(x, prec);

    mpfr_t sinval, cosval;

    mpfr_init2(sinval, prec);
    mpfr_init2(cosval, prec);

    for (unsigned i = 0; i < n; ++i)
    {
        mpfr_mul_ui(x, coeff, i, rnd);

        mpc_t * sig = new mpc_t[1];

        mpc_init2(sig[0], prec);

        mpfr_sin_cos(sinval, cosval, x, rnd);

        mpc_set_fr_fr(sig[0], sinval, cosval, rnd);

        sigvec.push_back(sig);
    }

    mpfr_clear(cosval);
    mpfr_clear(sinval);

    mpfr_clear(coeff);
    mpfr_clear(x);

    return sigvec;
}

vector<mpfr_t *> make_real_sinusoid(const mpfr_t & fs, const unsigned & n, const mpfr_t & f, const mpfr_prec_t & prec, const mpfr_rnd_t & rnd)
{
    vector<mpfr_t *> sigvec;
    sigvec.reserve(n);

    // calculate sin(2 * pi * f * i / fs);

    mpfr_t coeff;
    mpfr_init2(coeff, prec);

    mpfr_const_pi(coeff, rnd);
    mpfr_mul_ui(coeff, coeff, 2, rnd);
    mpfr_mul(coeff, coeff, f, rnd);
    mpfr_div(coeff, coeff, fs, rnd);

    //

    mpfr_t x;
    mpfr_init2(x, prec);

    for (unsigned i = 0; i < n; ++i)
    {
        mpfr_mul_ui(x, coeff, i, rnd);

        mpfr_t * sig = new mpfr_t[1];

        mpfr_init2(sig[0], prec);

        mpfr_sin(sig[0], x, rnd);

        sigvec.push_back(sig);
    }

    mpfr_clear(coeff);
    mpfr_clear(x);

    return sigvec;
}

int main()
{
    mpfr_t fs, f;

    mpfr_init2(fs, default_prec);
    mpfr_init2(f, default_prec);

    mpfr_set_d(fs, 2000.0, default_rnd);
    mpfr_set_d(f , 440.0 , default_rnd);

    vector<mpc_t *> sigvec = make_complex_sinusoid(fs, 100, f, default_prec, default_rnd);

    cout << "signal: ";

    {
        bool first = true;

        cout << '{';

        for (mpc_t * & ptr : sigvec)
        {
            if (first)
            {
                first = false;
            }
            else
            {
                cout << ", ";
            }

            char * str = mpc_get_str(10, 1000, *ptr, default_rnd);
            cout << str;
            mpc_free_str(str);
        }
    }

    cout << '}' << endl;

    for (mpc_t * & ptr : sigvec)
    {
        mpc_clear(*ptr);
        delete ptr;
    }

    mpfr_clear(fs);
    mpfr_clear(f);

    mpfr_free_cache();

    return 0;
}
