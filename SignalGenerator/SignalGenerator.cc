
////////////////////////
// SignalGenerator.cc //
////////////////////////

#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>

#if 0
const mpfr_rnd_t DEFAULT_MPFR_ROUNDINGMODE = MPFR_RNDN;
const mpc_rnd_t  DEFAULT_MPC_ROUNDINGMODE  = MPC_RNDNN;
#endif


void zero_signal(const double & x, double & y)
{
    // Name ........... : zero signal
    // Period ......... : n/a
    // DC offset ...... : 0
    // RMS ............ : 0
    // Range .......... : [ 0 ]

    (void)x; // unused
    y = 0;
}

void dc_signal(const double & x, double & y)
{
    // Name ........... : DC signal
    // Period ......... : n/a
    // DC offset ...... : sqrt(2) / 2
    // RMS ............ : sqrt(2) / 2
    // Range .......... : [ sqrt(2) / 2 ]

    (void)x; // unused
    y = sqrt(0.5);
}

void sine_signal(const double & x, double & y)
{
    // Name ........... : sine wave
    // Period ......... : 2*pi
    // DC offset ...... : 0
    // RMS ............ : sqrt(1 /2)
    // Range .......... : [ sqrt(2) / 2 ]

    y = sin(x);
}

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

using namespace std;

void get_info(const char * name, void f(const double & x, double & y), unsigned n)
{
    double min_value = +numeric_limits<double>::infinity();
    double max_value = -numeric_limits<double>::infinity();
    double sum1 = 0.0;
    double sum2 = 0.0;

    for (unsigned i = 0; i < n; ++i)
    {
        double x = static_cast<double>(i) / n * 2 * M_PI;
        double y;

        f(x, y);

        min_value = std::min(min_value, y);
        max_value = std::max(max_value, y);

        sum1 += y;
        sum2 += y * y;
    }

    double mean   = sum1 / n;
    double stddev = sqrt(abs(n * sum2 - sum1 * sum1)) / n; // The "abs" is there just to be sure we don't get to negative values due to roundoff.

    double power = sqrt(sum2 / n);

    cout << "================= " << name << endl;
    cout << endl;
    cout << "min value ....... : " << min_value       << endl;
    cout << "max value ....... : " << max_value       << endl;
    cout << "n ............... : " << n               << endl;
    cout << "sum1 ............ : " << sum1            << endl;
    cout << "sum2 ............ : " << sum2            << endl;
    cout << "mean value ...... : " << mean            << endl;
    cout << "stddev .......... : " << stddev          << endl;
    cout << "stddev^2 ........ : " << stddev * stddev << endl;
    cout << "power ........... : " << power           << endl;
    cout << "power^2 ......... : " << power * power   << endl;
    cout << endl;
}

int main()
{
    unsigned n = 10000000;

    get_info("sine_signal"     , sine_signal     , n);
    get_info("zero_signal"     , zero_signal     , n);
    get_info("dc_signal"       , dc_signal       , n);
    get_info("square_signal"   , square_signal   , n);
    get_info("triangle_signal" , triangle_signal , n);
    get_info("sawtooth_signal" , sawtooth_signal , n);
    get_info("gaussian_signal" , gaussian_signal , n);

    return 0;
}
