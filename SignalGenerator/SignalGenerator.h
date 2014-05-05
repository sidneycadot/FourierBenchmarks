
///////////////////////
// SignalGenerator.h //
///////////////////////

#ifndef SignalGenerator_h
#define SignalGenerator_h

#include <vector>
#include <complex>
#include <string>

#include <mpc.h>

#include "../MultiPrecisionUtils/MultiPrecisionUtils.h"

class GaussianNoiseSignal
{
    public:

        GaussianNoiseSignal(const std::string & seed_string, const mpfr_prec_t precision) :
            seed_string(seed_string)
        {
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
        }

        void operator () (mpc_t & rop, const std::vector<signed> & index)
        {
            // Use index to set up seed, so we're reproducible.

            std::string index_string = seed_string;

            for (unsigned i = 0; i < index.size(); ++i)
            {
                if (index_string.empty())
                {
                    // add separator.
                    index_string += ".";
                }

                index_string += std::to_string(index[i]);
            }

            gmp_randseed_string(randstate, index_string);

            mpfr_grandom(re, im, randstate, DEFAULT_MPFR_ROUNDINGMODE);

            mpc_set_fr_fr(rop, re, im, DEFAULT_MPFR_ROUNDINGMODE);
        }

    private:

        std::string seed_string;

        gmp_randstate_t randstate;

        mpfr_t re;
        mpfr_t im;
};

// NOTE: This version does NOT have reproducible values for an index.
class GaussianNoiseSignalFast
{
    public:

        // Note: reproducibility of randomly generated numbers is only ensured if precision is identical.

        GaussianNoiseSignalFast(const std::string & seed_string, const mpfr_prec_t precision)
        {
            // Prepare seed

            // Prepare "re" and "im" variables, used during random generation.

            mpfr_init2(re, precision);
            mpfr_init2(im, precision);

            gmp_randinit_default(randstate);

            gmp_randseed_string(randstate, seed_string);
        }

        ~GaussianNoiseSignalFast()
        {
            gmp_randclear(randstate);

            mpfr_clear(im);
            mpfr_clear(re);
        }

        void operator () (mpc_t & rop, const std::vector<signed> & index)
        {
            (void)index; // unused

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

// Routines to map an mpc_t to a complex or real value.

template <typename T>
void map_value(T & rop, const mpc_t & op)
{
    rop = mpfr_get_fp<T>(mpc_realref(op), DEFAULT_MPFR_ROUNDINGMODE);
}

template <typename T>
void map_value(std::complex<T> & rop, const mpc_t & op)
{
    rop = std::complex<T>(mpfr_get_fp<T>(mpc_realref(op), DEFAULT_MPFR_ROUNDINGMODE), mpfr_get_fp<T>(mpc_imagref(op), DEFAULT_MPFR_ROUNDINGMODE));
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
