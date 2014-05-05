
///////////////////////////
// MultiPrecisionUtils.h //
///////////////////////////

#ifndef MultiPrecisionUtils_h
#define MultiPrecisionUtils_h

#include <string>

#include <mpfr.h>
#include <mpc.h>

const mpfr_rnd_t DEFAULT_MPFR_ROUNDINGMODE = MPFR_RNDN;
const mpc_rnd_t  DEFAULT_MPC_ROUNDINGMODE  = MPC_RNDNN;

std::string to_string(const mpfr_t & x);

template <typename T>
inline T mpfr_to_native(const mpfr_t & op);

template <>
inline float mpfr_to_native<float>(const mpfr_t & op)
{
    return mpfr_get_flt(op, DEFAULT_MPFR_ROUNDINGMODE);
}

template <>
inline double mpfr_to_native<double>(const mpfr_t & op)
{
    return mpfr_get_d(op, DEFAULT_MPFR_ROUNDINGMODE);
}

template <>
inline long double mpfr_to_native<long double>(const mpfr_t & op)
{
    return mpfr_get_ld(op, DEFAULT_MPFR_ROUNDINGMODE);
}

#endif // MultiPrecisionUtils_h
