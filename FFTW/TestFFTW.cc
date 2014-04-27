
/////////////////
// TestFFTW.cc //
/////////////////

#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <chrono>

#include <fftw3.h>

using namespace std;

template <typename fp_t> struct FFTWX;

template <>
struct FFTWX<float>
{
    typedef float         real;
    typedef fftwf_complex complex;

    typedef fftwf_plan plan;

    static void free(void *p)
    {
        fftwf_free(p);
    }

    static complex * alloc_complex(size_t n)
    {
        return fftwf_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex * in, complex * out, int sign, unsigned flags)
    {
        return fftwf_plan_dft_1d(n0, in, out, sign, flags);
    }

    static void execute(const plan plan)
    {
        fftwf_execute(plan);
    }

    static void destroy_plan(plan plan)
    {
        fftwf_destroy_plan(plan);
    }

    static void cleanup(void)
    {
        fftwf_cleanup();
    }
};

template <>
struct FFTWX<double>
{
    typedef double       real;
    typedef fftw_complex complex;

    typedef fftw_plan plan;

    static void free(void *p)
    {
        fftw_free(p);
    }

    static complex * alloc_complex(size_t n)
    {
        return fftw_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex * in, complex * out, int sign, unsigned flags)
    {
        return fftw_plan_dft_1d(n0, in, out, sign, flags);
    }

    static void execute(const plan plan)
    {
        fftw_execute(plan);
    }

    static void destroy_plan(plan plan)
    {
        fftw_destroy_plan(plan);
    }

    static void cleanup(void)
    {
        fftw_cleanup();
    }
};

template <>
struct FFTWX<long double>
{
    typedef long double real;
    typedef fftwl_complex complex;

    typedef fftwl_plan plan;

    static void free(void *p)
    {
        fftwl_free(p);
    }

    static complex * alloc_complex(size_t n)
    {
        return fftwl_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex * in, complex * out, int sign, unsigned flags)
    {
        return fftwl_plan_dft_1d(n0, in, out, sign, flags);
    }

    static void execute(const plan plan)
    {
        fftwl_execute(plan);
    }

    static void destroy_plan(plan plan)
    {
        fftwl_destroy_plan(plan);
    }

    static void cleanup(void)
    {
        fftwl_cleanup();
    }
};

template <>
struct FFTWX<__float128>
{
    typedef __float128 real;
    typedef fftwq_complex complex;

    typedef fftwq_plan plan;

    static void free(void *p)
    {
        fftwq_free(p);
    }

    static complex * alloc_complex(size_t n)
    {
        return fftwq_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex * in, complex * out, int sign, unsigned flags)
    {
        return fftwq_plan_dft_1d(n0, in, out, sign, flags);
    }

    static void execute(const plan plan)
    {
        fftwq_execute(plan);
    }

    static void destroy_plan(plan plan)
    {
        fftwq_destroy_plan(plan);
    }

    static void cleanup(void)
    {
        fftwq_cleanup();
    }
};

template <typename fftwx>
void execute_test_c2c(const unsigned & n)
{
    const unsigned repeats = 5;

    vector<complex<typename fftwx::real>> v(n);

    typename fftwx::complex * x = fftwx::alloc_complex(n);
    assert(x != nullptr);

    typename fftwx::complex * y = fftwx::alloc_complex(n);
    assert(y != nullptr);

    typename fftwx::plan plan = fftwx::plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    assert(plan != nullptr);

    // Initialize the input

    // copy test data
    for (unsigned i = 0; i < n; ++i)
    {
        x[i][0] = real(v[i]);
        x[i][1] = imag(v[i]);
    }

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        auto t1 = chrono::high_resolution_clock::now();
        fftwx::execute(plan);
        auto t2 = chrono::high_resolution_clock::now();

        double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

        cout << "duration: " << duration << endl;
    }

    fftwx::destroy_plan(plan);

    fftwx::free(y);
    fftwx::free(x);

    fftwx::cleanup();
}

int main()
{
    const unsigned n = 65536;

    execute_test_c2c<FFTWX<float>>(n);
    execute_test_c2c<FFTWX<double>>(n);
    execute_test_c2c<FFTWX<long double>>(n);
    //execute_test_c2c<FFTWX<__float128>>(n);

    return EXIT_SUCCESS;
}
