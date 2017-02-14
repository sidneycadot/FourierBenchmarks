
////////////
// main.c //
////////////

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#include "fft_v0.h"
#include "fft_v1.h"
#include "czt.h"

typedef void fftfunc_t(complex double *, unsigned);

static void czt_fft_v1(complex double * z, unsigned size)
{
    czt_fft(z, size, fft_v1);
}

static void prepare_input(complex double * z, unsigned size)
{
    for (unsigned k = 0; k < size; ++k)
    {
        z[k] = cos(1 + k * k) + sin(1 + k * k * k) * I;
    }
}

static void demo(const char * fftname, fftfunc_t fftfunc, unsigned size)
{
    complex double z[size];

    prepare_input(z, size);

    fftfunc(z, size);

    for (unsigned k = 0; k < size; ++k)
    {
        printf("%-20s %10u %25.20f %25.20f\n", fftname, k, creal(z[k]), cimag(z[k]));
    }
    printf("\n");
}

static void time_fftfunc(const char * fftname, fftfunc_t fftfunc, unsigned size)
{
    complex double z[size];

    prepare_input(z, size);

    const unsigned num_warmup  = 9;
    const unsigned num_average = 16;

    const unsigned num_rep = num_warmup + num_average;

    clock_t timings[num_rep];

    for (unsigned rep = 0; rep < num_rep; ++rep)
    {
        clock_t t0 = clock();
        fftfunc(z, size);
        clock_t t1 = clock();

        timings[rep] = (t1 - t0);
    }

    clock_t total_time = 0;

    for (unsigned rep = num_warmup; rep < num_rep; ++rep)
    {
        total_time += timings[rep];
    }

    const double avg_time = (double)total_time / (double)(num_average * CLOCKS_PER_SEC);

    printf("%-20s %5u %20.9f\n", fftname, size, avg_time);
}

void time_fftfunc_loop(const char * fftname, fftfunc_t fftfunc, unsigned max_size)
{
    for (unsigned size = 1; size <= max_size; size *= 2)
    {
        time_fftfunc(fftname, fftfunc, size);
    }
    printf("\n");
}


int main()
{
    demo("fft_v0"    , fft_v0    , 5);
    demo("czt_fft_v1", czt_fft_v1, 5);

    time_fftfunc("czt_fft_v1", czt_fft_v1, 12345);

    //time_fftfunc_loop("fft_v0", fft_v0, 256);
    //time_fftfunc_loop("fft_v1", fft_v1, 16384);
    return 0;
}

