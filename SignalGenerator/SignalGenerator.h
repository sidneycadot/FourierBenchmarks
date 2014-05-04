
///////////////////////
// SignalGenerator.h //
///////////////////////

#ifndef SignalGenerator_h
#define SignalGenerator_h

#include <vector>
#include <complex>
#include <string>

#include <mpc.h>

const mpfr_rnd_t DEFAULT_MPFR_ROUNDINGMODE = MPFR_RNDN;
const mpc_rnd_t  DEFAULT_MPC_ROUNDINGMODE  = MPC_RNDNN;

class GaussianNoiseSignal
{
    public:

        GaussianNoiseSignal(const std::string & seed_string, const mpfr_prec_t precision)
        {
            // Prepare seed

            mpz_init(seed); // initialize and set to zero.
            mpz_init(seed2);

            for (unsigned i = 0; seed_string[i] != '\0'; ++i)
            {
                const unsigned c = static_cast<unsigned char>(seed_string[i]);

                mpz_mul_ui(seed, seed, 256);
                mpz_add_ui(seed, seed, c);
            }

            // Prepare "re" and "im" variables, used during random generation.

            mpfr_init2(re, precision);
            mpfr_init2(im, precision);

            gmp_randinit_default(randstate);
        }

        ~GaussianNoiseSignal()
        {
            gmp_randclear(randstate);

            mpfr_clear(im);
            mpfr_clear(re);

            mpz_clear(seed2);
            mpz_clear(seed);
        }

        void operator () (mpc_t & rop, const std::vector<signed> & index)
        {
            // Use index to set up seed, so we're reproducible.

            std::string index_string;
            for (unsigned i = 0; i < index.size(); ++i)
            {
                index_string += ".";
                index_string += std::to_string(index[i]);
            }

            mpz_set(seed2, seed);

            for (unsigned i = 0; index_string[i] != '\0'; ++i)
            {
                const unsigned c = static_cast<unsigned char>(index_string[i]);

                mpz_mul_ui(seed2, seed2, 256);
                mpz_add_ui(seed2, seed2, c);
            }

            gmp_randseed(randstate, seed2);

            mpfr_grandom(re, im, randstate, DEFAULT_MPFR_ROUNDINGMODE);

            mpc_set_fr_fr(rop, re, im, DEFAULT_MPFR_ROUNDINGMODE);
        }

    private:

        mpz_t seed;
        mpz_t seed2;

        gmp_randstate_t randstate;

        mpfr_t re;
        mpfr_t im;
};

// NOTE: This version does NOT have reproducible values for an index.
class GaussianNoiseSignalFast
{
    public:

        GaussianNoiseSignalFast(const std::string & seed_string, const mpfr_prec_t precision)
        {
            // Prepare seed

            mpz_t seed;

            mpz_init(seed); // initialize and set to zero.

            for (unsigned i = 0; seed_string[i] != '\0'; ++i)
            {
                const unsigned c = static_cast<unsigned char>(seed_string[i]);

                mpz_mul_ui(seed, seed, 256);
                mpz_add_ui(seed, seed, c);
            }

            // Prepare "re" and "im" variables, used during random generation.

            mpfr_init2(re, precision);
            mpfr_init2(im, precision);

            gmp_randinit_default(randstate);

            gmp_randseed(randstate, seed);

            mpz_clear(seed);
        }

        ~GaussianNoiseSignalFast()
        {
            gmp_randclear(randstate);

            mpfr_clear(im);
            mpfr_clear(re);
        }

        void operator () (mpc_t & rop, const std::vector<signed> & index)
        {
            (void)index;

            // Use index to set up seed, so we're reproducible.

            mpfr_grandom(re, im, randstate, DEFAULT_MPFR_ROUNDINGMODE);

            mpc_set_fr_fr(rop, re, im, DEFAULT_MPFR_ROUNDINGMODE);
        }

    private:

        gmp_randstate_t randstate;

        mpfr_t re;
        mpfr_t im;
};

class ZeroSignal
{
    public:

        void operator () (mpc_t & rop, const std::vector<signed> & index)
        {
            (void)index; // unused

            mpc_set_ui(rop, 0, DEFAULT_MPC_ROUNDINGMODE);
        }
};

template <typename T>
T mpfr_value(const mpfr_t & op);

template <>
float mpfr_value<float>(const mpfr_t & op)
{
    return mpfr_get_flt(op, DEFAULT_MPFR_ROUNDINGMODE);
}

template <>
double mpfr_value<double>(const mpfr_t & op)
{
    return mpfr_get_d(op, DEFAULT_MPFR_ROUNDINGMODE);
}

template <>
long double mpfr_value<long double>(const mpfr_t & op)
{
    return mpfr_get_ld(op, DEFAULT_MPFR_ROUNDINGMODE);
}

// Routines to map an mpc_t to a complex or real value.

template <typename T>
void map_value(T & rop, const mpc_t & op)
{
    rop = mpfr_value<T>(mpc_realref(op));
}

template <typename T>
void map_value(std::complex<T> & rop, const mpc_t & op)
{
    rop = std::complex<T>(mpfr_value<T>(mpc_realref(op)), mpfr_value<T>(mpc_imagref(op)));
}

template <typename Function, typename T>
void sample(Function & f, T * base, const std::vector<signed> & strides, const std::vector<unsigned> & dims, const mpfr_prec_t & precision)
{
    mpc_t value;

    mpc_init2(value, precision);

    std::vector<int> index(dims.size());

    long signed offset = 0;

    for (;;) // loop over the elements
    {
        // "index" and "offset" are now valid.

         // set value to f(index)

        f(value, index);

        map_value(base[offset], value);

        // Proceed to next item.

        int incdim = dims.size();

        while ((--incdim) >= 0)
        {
            ++index[incdim];
            offset += strides[incdim];

            if (static_cast<unsigned>(index[incdim]) != dims[incdim])
            {
                break;
            }

            // Proceed to earlier dimension

            index[incdim] = 0;
            offset -= strides[incdim] * dims[incdim];
        }

        if (incdim < 0)
        {
            // This was the last element
            break;
        }
    }

    mpc_clear(value);
}

#endif // SignalGenerator_h
