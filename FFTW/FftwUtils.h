
/////////////////
// FftwUtils.h //
/////////////////

#ifndef FftwUtil_h
#define FftwUtil_h

#include <fftw3.h>

template <typename fp_t> struct FFTW_Traits;

template <>
struct FFTW_Traits<float>
{
    typedef float         real_type;
    typedef fftwf_complex complex_type;

    typedef fftwf_plan plan;

    static void free(void *p)
    {
        fftwf_free(p);
    }

    static real_type * alloc_real(size_t n)
    {
        return fftwf_alloc_real(n);
    }

    static complex_type * alloc_complex(size_t n)
    {
        return fftwf_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex_type * in, complex_type * out, int sign, unsigned flags)
    {
        return fftwf_plan_dft_1d(n0, in, out, sign, flags);
    }

    static plan plan_dft_r2c_1d(int n, real_type * in, complex_type * out, unsigned flags)
    {
        return fftwf_plan_dft_r2c_1d(n, in, out, flags);
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
struct FFTW_Traits<double>
{
    typedef double       real_type;
    typedef fftw_complex complex_type;

    typedef fftw_plan plan;

    static void free(void *p)
    {
        fftw_free(p);
    }

    static real_type * alloc_real(size_t n)
    {
        return fftw_alloc_real(n);
    }

    static complex_type * alloc_complex(size_t n)
    {
        return fftw_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex_type * in, complex_type * out, int sign, unsigned flags)
    {
        return fftw_plan_dft_1d(n0, in, out, sign, flags);
    }

    static plan plan_dft_r2c_1d(int n, real_type * in, complex_type * out, unsigned flags)
    {
        return fftw_plan_dft_r2c_1d(n, in, out, flags);
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
struct FFTW_Traits<long double>
{
    typedef long double real_type;
    typedef fftwl_complex complex_type;

    typedef fftwl_plan plan;

    static void free(void *p)
    {
        fftwl_free(p);
    }

    static real_type * alloc_real(size_t n)
    {
        return fftwl_alloc_real(n);
    }

    static complex_type * alloc_complex(size_t n)
    {
        return fftwl_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex_type * in, complex_type * out, int sign, unsigned flags)
    {
        return fftwl_plan_dft_1d(n0, in, out, sign, flags);
    }

    static plan plan_dft_r2c_1d(int n, real_type * in, complex_type * out, unsigned flags)
    {
        return fftwl_plan_dft_r2c_1d(n, in, out, flags);
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
struct FFTW_Traits<__float128>
{
    typedef __float128    real_type;
    typedef fftwq_complex complex_type;

    typedef fftwq_plan plan;

    static void free(void *p)
    {
        fftwq_free(p);
    }

    static real_type * alloc_real(size_t n)
    {
        return fftwq_alloc_real(n);
    }

    static complex_type * alloc_complex(size_t n)
    {
        return fftwq_alloc_complex(n);
    }

    static plan plan_dft_1d(int n0, complex_type * in, complex_type * out, int sign, unsigned flags)
    {
        return fftwq_plan_dft_1d(n0, in, out, sign, flags);
    }

    static plan plan_dft_r2c_1d(int n, real_type * in, complex_type * out, unsigned flags)
    {
        return fftwq_plan_dft_r2c_1d(n, in, out, flags);
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

#endif // FftwUtil_h
