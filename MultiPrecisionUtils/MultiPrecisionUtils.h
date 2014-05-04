
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

#endif // MultiPrecisionUtils_h
