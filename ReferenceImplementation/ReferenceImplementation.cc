
////////////////////////////////
// ReferenceImplementation.cc //
////////////////////////////////

#include "ReferenceImplementation.h"

#include <cassert>

void pow2_fft(const FourierTransformDirection direction, mpc_t * z, const unsigned n, const unsigned stride, const mpfr_prec_t precision)
{
    if (n == 1)
    {
        return;
    }

    assert (n % 2 == 0);

    // Do sub-FFTs on even / odd entries.

    pow2_fft(direction, z         , n / 2, 2 * stride, precision);
    pow2_fft(direction, z + stride, n / 2, 2 * stride, precision);

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
            mpfr_mul_si(turn, turn, -2 * i, DEFAULT_MPFR_ROUNDINGMODE);
        }
        else
        {
            mpfr_mul_ui(turn, turn, +2 * i, DEFAULT_MPFR_ROUNDINGMODE);
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

void czt(const mpc_t * x, unsigned n, mpc_t * y, unsigned m, const mpc_t & w, const mpc_t & a, const mpfr_prec_t precision)
{
    // Make sure we have room for (m + n - 1) entries
    // Determine next-biggest power-of-two

    unsigned fft_size = 1;
    while (fft_size < (m + n - 1))
    {
        fft_size += fft_size;
    }

    mpc_t temp;
    mpc_init2(temp, precision);

    // xx = x * a ** -np.arange(n) * w ** np.arange(0, n) ** 2 (padded by zeroes)

    mpc_t * xx = new mpc_t[fft_size];

    for (unsigned i = 0; i < fft_size; ++i)
    {
        mpc_init2(xx[i], precision);
    }

    // xx = x * a ** -np.arange(n) * w ** wExponents1

    for (unsigned i = 0; i < n; ++i)
    {
        // xx[i] = x[i] * pow(w, (i * i) / 2) / pow(a, i)

        mpc_set_ui(temp, (i * i), DEFAULT_MPC_ROUNDINGMODE);
        mpc_div_ui(temp, temp, 2, DEFAULT_MPC_ROUNDINGMODE);

        mpc_pow(temp, w, temp, DEFAULT_MPC_ROUNDINGMODE);

        mpc_mul(xx[i], x[i], temp, DEFAULT_MPC_ROUNDINGMODE); // xx[i] = x[i] 

        mpc_pow_ui(temp, a, i, DEFAULT_MPC_ROUNDINGMODE);

        mpc_div(xx[i], xx[i], temp, DEFAULT_MPC_ROUNDINGMODE);
    }

    for (unsigned i = n; i < fft_size; ++i) // pad with zeroes
    {
        mpc_set_ui(xx[i], 0, DEFAULT_MPC_ROUNDINGMODE);
    }

    // Do FFT of xx

    pow2_fft(FourierTransformDirection::Forward, xx, fft_size, 1, precision);

    // Allocate & initialize ww

    mpc_t * ww = new mpc_t[fft_size];

    for (unsigned i = 0; i < fft_size; ++i)
    {
        mpc_init2(ww[i], precision);
    }

    for (unsigned i = 0; i < m + n - 1; ++i)
    {
        // ww[i] = pow(w, -(j * j) / 2)  with j = i - (n - 1)

        int j = i - (n - 1);

        mpc_set_si(temp, -(j * j), DEFAULT_MPC_ROUNDINGMODE);
        mpc_div_ui(temp, temp, 2, DEFAULT_MPC_ROUNDINGMODE);

        mpc_pow(ww[i], w, temp, DEFAULT_MPC_ROUNDINGMODE);
    }

    for (unsigned i = m + n - 1; i < fft_size; ++i) // pad with zeroes
    {
        mpc_set_ui(ww[i], 0, DEFAULT_MPC_ROUNDINGMODE);
    }

    // Do FFT of ww

    pow2_fft(FourierTransformDirection::Forward, ww, fft_size, 1, precision);

    // Do convolution: xx = (xx * ww), element-wise

    for (unsigned i = 0; i < fft_size; ++i)
    {
        mpc_mul(xx[i], xx[i], ww[i], DEFAULT_MPC_ROUNDINGMODE);
    }

    // Do inverse FFT of (xx * ww), put result in xx.

    pow2_fft(FourierTransformDirection::Inverse, xx, fft_size, 1, precision);

    // The last 'm' elements are the output (after multiplying)
    // Put them in 'y'.

    for (unsigned i = 0; i < m; ++i)
    {
        mpc_set_ui(temp, (i * i), DEFAULT_MPC_ROUNDINGMODE);
        mpc_div_ui(temp, temp, 2, DEFAULT_MPC_ROUNDINGMODE);

        mpc_pow(temp, w, temp, DEFAULT_MPC_ROUNDINGMODE);

        mpc_mul(y[i], xx[i + n - 1], temp, DEFAULT_MPC_ROUNDINGMODE);
    }

    // Deallocate ww, xx, and temp

    for (unsigned i = 0; i < fft_size; ++i)
    {
        mpc_clear(ww[i]);
    }

    delete [] ww;

    for (unsigned i = 0; i < fft_size; ++i)
    {
        mpc_clear(xx[i]);
    }

    delete [] xx;

    mpc_clear(temp);
}

void generic_fft(const FourierTransformDirection direction, mpc_t * z, const unsigned n, const mpfr_prec_t precision)
{
    const unsigned & m = n;

    // Init w
    mpc_t w;
    mpc_init2(w, precision);

    {
        // w = np.exp(-2j * np.pi / m)

        mpfr_t turn;
        mpfr_t sin_turn;
        mpfr_t cos_turn;

        mpfr_init2(turn    , precision);
        mpfr_init2(sin_turn, precision);
        mpfr_init2(cos_turn, precision);

        mpfr_const_pi(turn, DEFAULT_MPFR_ROUNDINGMODE);

        if (direction == FourierTransformDirection::Forward)
        {
            mpfr_mul_si(turn, turn, -2, DEFAULT_MPFR_ROUNDINGMODE);
        }
        else
        {
            mpfr_mul_si(turn, turn, +2, DEFAULT_MPFR_ROUNDINGMODE);
        }

        mpfr_div_ui(turn, turn, m, DEFAULT_MPFR_ROUNDINGMODE);

        // w = exp(turn * i)

        mpfr_sin_cos(sin_turn, cos_turn, turn, DEFAULT_MPFR_ROUNDINGMODE);

        mpc_set_fr_fr(w, cos_turn, sin_turn, DEFAULT_MPC_ROUNDINGMODE);

        mpfr_clear(cos_turn);
        mpfr_clear(sin_turn);
        mpfr_clear(turn);
    }

    // Init a

    mpc_t a;
    mpc_init2(a, precision);

    mpc_set_ui(a, 1, DEFAULT_MPC_ROUNDINGMODE);

    // Do FFT as CZT

    czt(z, n, z, m, w, a, precision);

    if (direction == FourierTransformDirection::Inverse)
    {
        for (unsigned i = 0; i < m; ++i)
        {
            mpc_div_ui(z[i], z[i], m, DEFAULT_MPC_ROUNDINGMODE);
        }
    }

    // Remove old stuff

    mpc_clear(a);
    mpc_clear(w);
}
