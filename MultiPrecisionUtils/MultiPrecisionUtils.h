
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

// get native floating-point type from mpfr

template <typename T>
inline T mpfr_get_fp(const mpfr_t & op, const mpfr_rnd_t rnd);

template <>
inline float mpfr_get_fp<float>(const mpfr_t & op, const mpfr_rnd_t rnd)
{
    return mpfr_get_flt(op, rnd);
}

template <>
inline double mpfr_get_fp<double>(const mpfr_t & op, const mpfr_rnd_t rnd)
{
    return mpfr_get_d(op, rnd);
}

template <>
inline long double mpfr_get_fp<long double>(const mpfr_t & op, const mpfr_rnd_t rnd)
{
    return mpfr_get_ld(op, rnd);
}

// assign mpc from native floating-point type (real part only)

inline void mpc_set_fp(mpc_t & rop, const float & re, const mpc_rnd_t rnd)
{
    // Note: there is no set with explicit 'float' arguments.
    mpc_set_d(rop, re, rnd);
}

inline void mpc_set_fp(mpc_t & rop, const double & re, const mpc_rnd_t rnd)
{
    mpc_set_d(rop, re, rnd);
}

inline void mpc_set_fp(mpc_t & rop, const long double & re, const mpc_rnd_t rnd)
{
    mpc_set_ld(rop, re, rnd);
}

// assign mpc from native floating-point type (real and imaginary parts)

inline void mpc_set_fp_fp(mpc_t & rop, const float & re, const float & im, const mpc_rnd_t rnd)
{
    // Note: there is no set with explicit 'float' arguments.
    mpc_set_d_d(rop, re, im, rnd);
}

inline void mpc_set_fp_fp(mpc_t & rop, const double & re, const double & im, const mpc_rnd_t rnd)
{
    mpc_set_d_d(rop, re, im, rnd);
}

inline void mpc_set_fp_fp(mpc_t & rop, const long double & re, const long double & im, const mpc_rnd_t rnd)
{
    mpc_set_ld_ld(rop, re, im, rnd);
}

#endif // MultiPrecisionUtils_h
