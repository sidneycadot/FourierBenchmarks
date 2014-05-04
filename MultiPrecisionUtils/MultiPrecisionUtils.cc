
////////////////////////////
// MultiPrecisionUtils.cc //
////////////////////////////

#include <cassert>

#include "MultiPrecisionUtils.h"

std::string to_string(const mpfr_t & x)
{
    char * str;
    int mpfr_asprintf_result = mpfr_asprintf(&str, "%.6Re", x);
    assert(mpfr_asprintf_result > 0);

    std::string r = str;

    mpfr_free_str(str);

    return r;
}
