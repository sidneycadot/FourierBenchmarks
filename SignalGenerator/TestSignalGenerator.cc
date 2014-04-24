
////////////////////////////
// TestSignalGenerator.cc //
////////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cassert>
#include <string>

#include <mpfr.h>

const mpfr_rnd_t DEFAULT_MPFR_ROUNDINGMODE = MPFR_RNDN;

#if 0
const mpc_rnd_t  DEFAULT_MPC_ROUNDINGMODE  = MPC_RNDNN;
#endif

std::string to_string(const mpfr_t & x)
{
    char * str;
    int mpfr_asprintf_result = mpfr_asprintf(&str, "%.6Re", x);
    assert(mpfr_asprintf_result > 0);

    std::string r = str;

    mpfr_free_str(str);

    return r;
}

void zero_signal(const mpfr_t & x, mpfr_t & y)
{
    // Name ........... : zero signal
    // Period ......... : n/a
    // DC offset ...... : 0
    // RMS ............ : 0
    // Range .......... : [ 0 ]

    (void)x; // unused

    mpfr_set_zero(y, +1); // y = +0
}


void dc_signal(const mpfr_t & x, mpfr_t & y)
{
    // Name ........... : DC signal
    // Period ......... : n/a
    // DC offset ...... : sqrt(2) / 2
    // RMS ............ : sqrt(2) / 2
    // Range .......... : [ sqrt(2) / 2 ]

    (void)x; // unused

    // y = sqrt(0.5);

    mpfr_set_ui(y, 1, DEFAULT_MPFR_ROUNDINGMODE);
    mpfr_div_ui(y, y, 2, DEFAULT_MPFR_ROUNDINGMODE);
    mpfr_sqrt(y, y, DEFAULT_MPFR_ROUNDINGMODE);
}

void sine_signal(const mpfr_t & x, mpfr_t & y)
{
    // Name ........... : sine wave
    // Period ......... : 2*pi
    // DC offset ...... : 0
    // RMS ............ : sqrt(1 /2)
    // Range .......... : [ sqrt(2) / 2 ]

    mpfr_sin(y, x, DEFAULT_MPFR_ROUNDINGMODE);
}

#if 0
void square_signal(const double & x, double & y)
{
    // Name ........... : square wave
    // Period ......... : 2*pi
    // DC offset ...... : 0
    // RMS ............ : sqrt(1 / 2)
    // Range .......... : [ -sqrt(1 / 2), +sqrt(1 / 2) ]

    double xnorm = x / (2 * M_PI);

    double frac = xnorm - floor(xnorm);

    y = sqrt(0.5) * ((frac < 0.5) ? +1 : -1);
}

void triangle_signal(const double & x, double & y)
{
    // Name ........... : triangle wave
    // Period ......... : 2*pi
    // DC offset ...... : 0
    // RMS ............ : sqrt(2)/2
    // Range .......... : [ -sqrt(3 / 2), +sqrt(3 / 2) ]

    double xnorm = x / (2 * M_PI);

    double frac = xnorm - floor(xnorm);

    y = sqrt(6) * (frac - abs(frac - 0.25) + abs(frac - 0.75) - 0.5);
}

void sawtooth_signal(const double & x, double & y)
{
    // Name ........... : sawtooth wave
    // Period ......... : 2*pi
    // DC offset ...... : 0
    // RMS ............ : sqrt(2)/2
    // Range .......... : [ -C, +C ]

    double xnorm = x / (2 * M_PI);

    double frac = (xnorm - floor(xnorm));

    y = sqrt(6) * (frac - 0.5);
}

void gaussian_signal(const double & x, double & y)
{
    // Name ........... : gaussian noise
    // Period ......... : n/a
    // DC offset ...... : 0
    // RMS ............ : sqrt(2)/2
    // Range .......... : <-inf, +inf>

    (void)x;

    double U = static_cast<double>(rand()) / RAND_MAX;
    double V = static_cast<double>(rand()) / RAND_MAX;

    double R1 = sqrt(-2 * log(U)) * cos(2 * M_PI * V);
    //double R2 = sqrt(-2 * log(U)) * sin(2 * M_PI * V);

    y = sqrt(0.5) * R1;
}
#endif

using namespace std;

void get_info(const char * name, void f(const mpfr_t & x, mpfr_t & y), unsigned n, mpfr_prec_t precision)
{
    mpfr_t min_value;
    mpfr_t max_value;
    mpfr_t sum1;
    mpfr_t sum2;
    mpfr_t x;
    mpfr_t y;
    mpfr_t mean;
    mpfr_t stddevSquared;
    mpfr_t stddev;
    mpfr_t powerSquared;
    mpfr_t power;

    mpfr_init2(min_value     , precision);
    mpfr_init2(max_value     , precision);
    mpfr_init2(sum1          , precision);
    mpfr_init2(sum2          , precision);
    mpfr_init2(x             , precision);
    mpfr_init2(y             , precision);
    mpfr_init2(mean          , precision);
    mpfr_init2(stddevSquared , precision);
    mpfr_init2(stddev        , precision);
    mpfr_init2(powerSquared  , precision);
    mpfr_init2(power         , precision);

    mpfr_set_inf(min_value, +1); // +infinity
    mpfr_set_inf(max_value, -1); // -infinity

    mpfr_set_zero(sum1, +1); // +zero
    mpfr_set_zero(sum2, +1); // +zero

    for (unsigned i = 0; i < n; ++i)
    {
        // x = (double)i / n * 2 * M_PI;

        mpfr_const_pi(x, DEFAULT_MPFR_ROUNDINGMODE);
        mpfr_mul_ui(x, x, 2 * i, DEFAULT_MPFR_ROUNDINGMODE);
        mpfr_div_ui(x, x, n, DEFAULT_MPFR_ROUNDINGMODE);

        f(x, y);

        mpfr_min(min_value, min_value, y, DEFAULT_MPFR_ROUNDINGMODE);
        mpfr_max(max_value, max_value, y, DEFAULT_MPFR_ROUNDINGMODE);

        // sum1 += y

        mpfr_add(sum1, sum1, y, DEFAULT_MPFR_ROUNDINGMODE);

        // sum2 += y * y

        mpfr_sqr(y, y, DEFAULT_MPFR_ROUNDINGMODE);
        mpfr_add(sum2, sum2, y, DEFAULT_MPFR_ROUNDINGMODE);
    }

    // double mean   = sum1 / n;

    mpfr_div_ui(mean, sum1, n, DEFAULT_MPFR_ROUNDINGMODE);

    //double stddev = sqrt(abs(n * sum2 - sum1 * sum1)) / n; // The "abs" is there just to be sure we don't get to negative values due to roundoff.
    // x = sum2/n
    // y = mean*mean
    // stddevSquared = x - y

    mpfr_div_ui(x, sum2, n, DEFAULT_MPFR_ROUNDINGMODE);
    mpfr_sqr(y, mean, DEFAULT_MPFR_ROUNDINGMODE);

    mpfr_sub(stddevSquared, x, y, DEFAULT_MPFR_ROUNDINGMODE);

    // Due to roundoff, this could be just below 0.
    // Take the abs() to prevent this problem.
    mpfr_abs(stddevSquared, stddevSquared, DEFAULT_MPFR_ROUNDINGMODE);

    mpfr_sqrt(stddev, stddevSquared, DEFAULT_MPFR_ROUNDINGMODE);

    // double power = sqrt(sum2 / n);

    mpfr_div_ui(powerSquared, sum2, n, DEFAULT_MPFR_ROUNDINGMODE);

    mpfr_sqrt(power, powerSquared, DEFAULT_MPFR_ROUNDINGMODE);

    cout << "================= " << name << endl;
    cout << endl;

    cout << "min value ....... : " << to_string(min_value)       << endl;
    cout << "max value ....... : " << to_string(max_value)       << endl;
    cout << "n ............... : " << to_string(n)               << endl;
    cout << "sum1 ............ : " << to_string(sum1)            << endl;
    cout << "sum2 ............ : " << to_string(sum2)            << endl;
    cout << "mean value ...... : " << to_string(mean)            << endl;
    cout << "stddev .......... : " << to_string(stddev)          << endl;
    cout << "stddev^2 ........ : " << to_string(stddevSquared)   << endl;
    cout << "power ........... : " << to_string(power)           << endl;
    cout << "power^2 ......... : " << to_string(powerSquared)    << endl;

    cout << endl;

    mpfr_clear( min_value     );
    mpfr_clear( max_value     );
    mpfr_clear( sum1          );
    mpfr_clear( sum2          );
    mpfr_clear( x             );
    mpfr_clear( y             );
    mpfr_clear( mean          );
    mpfr_clear( stddevSquared );
    mpfr_clear( stddev        );
    mpfr_clear( powerSquared  );
    mpfr_clear( power         );
}

int main()
{
    unsigned n = 100000;

    mpfr_prec_t precision = 256;

    get_info("zero_signal"     , zero_signal     , n, precision);
    get_info("dc_signal"       , dc_signal       , n, precision);
    get_info("sine_signal"     , sine_signal     , n, precision);
    //get_info("square_signal"   , square_signal   , n);
    //get_info("triangle_signal" , triangle_signal , n);
    //get_info("sawtooth_signal" , sawtooth_signal , n);
    //get_info("gaussian_signal" , gaussian_signal , n);

    // Free cache of MPFR, used for PI.
    mpfr_free_cache();

    // Free the GMP random state.
    //gmp_randclear(rnd_state);

    return 0;
}
