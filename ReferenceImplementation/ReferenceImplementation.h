
///////////////////////////////
// ReferenceImplementation.h //
///////////////////////////////

#ifndef ReferenceImplementation_h
#define ReferenceImplementation_h

#include <mpc.h>

const mpfr_rnd_t DEFAULT_MPFR_ROUNDINGMODE = MPFR_RNDN;
const mpc_rnd_t  DEFAULT_MPC_ROUNDINGMODE  = MPC_RNDNN;

enum class FourierTransformDirection
{
    Forward,
    Inverse
};

void pow2_fft(const FourierTransformDirection direction, mpc_t * z, const unsigned n, const unsigned stride, const mpfr_prec_t precision);

void czt(const mpc_t * x, unsigned n, mpc_t * y, unsigned m, const mpc_t & w, const mpc_t & a, const mpfr_prec_t precision);

void generic_fft(const FourierTransformDirection direction, mpc_t * z, const unsigned n, const mpfr_prec_t precision);

#endif // ReferenceImplementation_h
