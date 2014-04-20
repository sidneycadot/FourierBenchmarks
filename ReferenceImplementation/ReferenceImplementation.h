
///////////////////////////////
// ReferenceImplementation.h //
///////////////////////////////

#ifndef ReferenceImplementation_h
#define ReferenceImplementation_h

#include <cassert>
#include <mpc.h>

const mpfr_rnd_t DEFAULT_MPFR_ROUNDINGMODE = MPFR_RNDN;
const mpc_rnd_t  DEFAULT_MPC_ROUNDINGMODE  = MPC_RNDNN;

enum class FourierTransformDirection
{
    Forward,
    Inverse
};

template <FourierTransformDirection direction>
void pow2_recursive_fft(mpc_t * z, const unsigned n, const unsigned stride, const mpfr_prec_t precision)
{
    if (n == 1)
    {
        return;
    }

    assert (n % 2 == 0);

    // Do sub-FFTs on even / odd entries.

    pow2_recursive_fft<direction>(z         , n / 2, 2 * stride, precision);
    pow2_recursive_fft<direction>(z + stride, n / 2, 2 * stride, precision);

    // Now, we need to combine the values of the sub-FFTs

    // First, allocate temp[] working space.

    mpc_t * temp = new mpc_t[n];

    for (unsigned i = 0; i < n; ++i)
    {
        mpc_init2(temp[i], precision);
    }

    // Calculate temp[] values as sum/difference of sub-FFT values.

    mpfr_t turn;
    mpfr_t sin_turn;
    mpfr_t cos_turn;
    mpc_t  coeff;

    mpfr_init2(turn    , precision);
    mpfr_init2(sin_turn, precision);
    mpfr_init2(cos_turn, precision);
    mpc_init2 (coeff   , precision);

    for (unsigned i = 0; i < n / 2; ++i)
    {
        // turn = -2.0 * M_PI * i / n [forward transform]
        // turn = +2.0 * M_PI * i / n [inverse transform]

        mpfr_const_pi(turn, DEFAULT_MPFR_ROUNDINGMODE);

        if (direction == FourierTransformDirection::Forward)
        {
            mpfr_mul_si(turn, turn, -2 * i , DEFAULT_MPFR_ROUNDINGMODE);
        }
        else
        {
            mpfr_mul_ui(turn, turn, +2 * i , DEFAULT_MPFR_ROUNDINGMODE);
        }

        mpfr_div_ui(turn, turn, n, DEFAULT_MPFR_ROUNDINGMODE);

        // coeff = exp(turn * i)

        mpfr_sin_cos(sin_turn, cos_turn, turn, DEFAULT_MPFR_ROUNDINGMODE);

        mpc_set_fr_fr(coeff, cos_turn, sin_turn, DEFAULT_MPC_ROUNDINGMODE);

        // Forward transform:
        //
        // temp[2 * i + 0] = z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)];
        // temp[2 * i + 1] = z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)];
        //
        // Inverse transform:
        //
        // temp[2 * i + 0] = 0.5 * (z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)]);
        // temp[2 * i + 1] = 0.5 * (z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)]);

        mpc_mul(coeff, coeff, z[stride * (2 * i + 1)], DEFAULT_MPC_ROUNDINGMODE);

        mpc_add(temp[2 * i + 0], z[stride * 2 * i], coeff, DEFAULT_MPC_ROUNDINGMODE);
        mpc_sub(temp[2 * i + 1], z[stride * 2 * i], coeff, DEFAULT_MPC_ROUNDINGMODE);

        if (direction == FourierTransformDirection::Inverse)
        {
            mpc_div_ui(temp[2 * i + 0], temp[2 * i + 0], 2, DEFAULT_MPC_ROUNDINGMODE);
            mpc_div_ui(temp[2 * i + 1], temp[2 * i + 1], 2, DEFAULT_MPC_ROUNDINGMODE);
        }
    }

    mpc_clear (coeff);
    mpfr_clear(cos_turn);
    mpfr_clear(sin_turn);
    mpfr_clear(turn);

    // Rearrange the temp[] entries into the z[] array.

    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned j = (i / 2) + (i % 2) * (n / 2);

        assert(j < n);

        mpc_swap(z[j * stride], temp[i]);
    }

    // Deallocate the temp[] array

    for (unsigned i = 0; i < n; ++i)
    {
        mpc_clear(temp[i]);
    }

    delete [] temp;

    // All done.
}

#if 0
void czt(const mpz_t * x, unsigned n, mpz_t * y, unsigned m, const mpc_t & w, const mpc_t a)
{
    //wExponents1 =   np.arange(    0   , n) ** 2
    //wExponents2 = - np.arange(-(n - 1), m) ** 2
    //wExponents3 =   np.arange(    0   , m) ** 2

    xx = x * a ** -np.arange(n) * w ** wExponents1

    // Determine next-biggest FFT of power-of-two

    unsigned nfft = 1;
    while (nfft < (m + n - 1))
    {
        nfft += nfft;
    }

    fxx = np.fft.fft(xx, nfft)

    fw = np.fft.fft(w ** wExponents2, nfft)

    fyy = fxx * fw

    yy = np.fft.ifft(fyy, nfft)

    // select output

    //yy = yy[n - 1 : m + n - 1]

    //y = yy * w ** wExponents3
}
#endif

#endif // ReferenceImplementation_h
