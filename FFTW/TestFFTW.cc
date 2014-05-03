
/////////////////
// TestFFTW.cc //
/////////////////

#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cassert>
#include <stdexcept>
#include <chrono>
#include <array>
#include <algorithm>

#include <fftw3.h>

using namespace std;

template <typename fp_t> struct FFTW_Traits;

template <typename T>
struct TypeName {};

template <>
struct TypeName<float>
{
    static std::string name() { return "float"; }
};

template <>
struct TypeName<double>
{
    static std::string name() { return "double"; }
};

template <>
struct TypeName<long double>
{
    static std::string name() { return "long_double"; }
};

string flags_to_string(unsigned flags)
{
    string description; // empty string means: not found yet.

    for (unsigned effort = 0; effort <= 3; ++effort)
    {
        for (unsigned allow_destroy = 0; allow_destroy <= 1; ++allow_destroy)
        {
            for (unsigned assume_aligned = 0; assume_aligned <= 1; ++assume_aligned)
            {
                unsigned test_flags = 0;

                if (effort == 0) test_flags |= FFTW_ESTIMATE;
                if (effort == 1) test_flags |= FFTW_MEASURE;
                if (effort == 2) test_flags |= FFTW_PATIENT;
                if (effort == 3) test_flags |= FFTW_EXHAUSTIVE;

                if (allow_destroy == 0) test_flags |= FFTW_PRESERVE_INPUT;
                if (allow_destroy == 1) test_flags |= FFTW_DESTROY_INPUT;

                if (assume_aligned == 0) test_flags |= FFTW_UNALIGNED; // no assumption that data is well aligned for SIMD.
                if (assume_aligned == 1) test_flags |= 0;

                if (test_flags == flags)
                {
                    assert(description.empty());

                    if (effort == 0) description = "E 0";
                    if (effort == 1) description = "E 1";
                    if (effort == 2) description = "E 2";
                    if (effort == 3) description = "E 3";

                    if (allow_destroy == 0) description += " D 0";
                    if (allow_destroy == 1) description += " D 1";

                    if (assume_aligned == 0) description += " A 0";
                    if (assume_aligned == 1) description += " A 1";
                }
            }
        }
    }

    assert(!description.empty());
    return description;
}

enum CacheTemperature {
    CacheCold = 0,
    CacheWarm = 1
};

enum SignalType {
    SignalZero,
    SignalDC,
    SignalNoise
};

string cachetemp_to_string(const CacheTemperature & ct)
{
    if (ct == CacheCold)
    {
        return "0";
    }

    if (ct == CacheWarm)
    {
        return "1";
    }

    assert(false);
}

string sigtype_to_string(const SignalType & st)
{
    if (st == SignalZero)
    {
        return "ZERO";
    }

    if (st == SignalDC)
    {
        return "DC";
    }

    if (st == SignalNoise)
    {
        return "NOISE";
    }

    assert(false);
}

// This is taken from the documentation of FFTW 3.3.4.
//
// 4.3.2 Planner Flags
//
// All of the planner routines in FFTW accept an integer flags argument, which is a bitwise OR (‘|’)
// of zero or more of the flag constants defined below. These flags control the rigor (and time) of
// the planning process, and can also impose (or lift) restrictions on the type of transform
// algorithm that is employed.
// 
// Important: the planner overwrites the input array during planning unless a saved plan (see Wisdom)
//            is available for that problem, so you should initialize your input data after creating
//            the plan. The only exceptions to this are the FFTW_ESTIMATE and FFTW_WISDOM_ONLY flags,
//            as mentioned below.
//
// In all cases, if wisdom is available for the given problem that was created with equal-or-greater
// planning rigor, then the more rigorous wisdom is used. For example, in FFTW_ESTIMATE mode any
// available wisdom is used, whereas in FFTW_PATIENT mode only wisdom created in patient or exhaustive
// mode can be used. See Words of Wisdom-Saving Plans.
//
// Planning-rigor flags
// 
//     FFTW_ESTIMATE    specifies that, instead of actual measurements of different algorithms, a simple
//                      heuristic is used to pick a (probably sub-optimal) plan quickly. With this flag,
//                      the input/output arrays are not overwritten during planning.
//
//     FFTW_MEASURE     tells FFTW to find an optimized plan by actually computing several FFTs and measuring
//                      their execution time. Depending on your machine, this can take some time (often a few
//                      seconds). FFTW_MEASURE is the default planning option.
//
//     FFTW_PATIENT     is like FFTW_MEASURE, but considers a wider range of algorithms and often produces a
//                      “more optimal” plan (especially for large transforms), but at the expense of several
//                      times longer planning time (especially for large transforms).
//
//     FFTW_EXHAUSTIVE  is like FFTW_PATIENT, but considers an even wider range of algorithms, including many
//                      that we think are unlikely to be fast, to produce the most optimal plan but with a
//                      substantially increased planning time.
//
//     FFTW_WISDOM_ONLY is a special planning mode in which the plan is only created if wisdom is available
//                      for the given problem, and otherwise a NULL plan is returned. This can be combined
//                      with other flags, e.g. ‘FFTW_WISDOM_ONLY | FFTW_PATIENT’ creates a plan only if
//                      wisdom is available that was created in FFTW_PATIENT or FFTW_EXHAUSTIVE mode.
//
//                      The FFTW_WISDOM_ONLY flag is intended for users who need to detect whether wisdom is
//                      available; for example, if wisdom is not available one may wish to allocate new arrays
//                      for planning so that user data is not overwritten. 
//
// Algorithm-restriction flags
//
//     FFTW_DESTROY_INPUT  specifies that an out-of-place transform is allowed to overwrite its input array
//                         with arbitrary data; this can sometimes allow more efficient algorithms to be employed.
//     FFTW_PRESERVE_INPUT specifies that an out-of-place transform must not change its input array.
//                         This is ordinarily the default, except for c2r and hc2r (i.e. complex-to-real) transforms
//                         for which FFTW_DESTROY_INPUT is the default. In the latter cases, passing
//                         FFTW_PRESERVE_INPUT will attempt to use algorithms that do not destroy the input, at the
//                         expense of worse performance; for multi-dimensional c2r transforms, however,
//                         no input-preserving algorithms are implemented and the planner will return NULL if one is requested.
//     FFTW_UNALIGNED      specifies that the algorithm may not impose any unusual alignment requirements on the input/output
//                         arrays (i.e. no SIMD may be used). This flag is normally not necessary, since the planner
//                         automatically detects misaligned arrays. The only use for this flag is if you want to use the new-array
//                         execute interface to execute a given plan on a different array that may not be aligned like the original.
//                         (Using fftw_malloc makes this flag unnecessary even then. You can also use fftw_alignment_of to detect
//                         whether two arrays are equivalently aligned.) 
//
// Limiting planning time
//
//      extern void fftw_set_timelimit(double seconds);
//
// This function instructs FFTW to spend at most seconds seconds (approximately) in the planner. If seconds == FFTW_NO_TIMELIMIT
// (the default value, which is negative), then planning time is unbounded. Otherwise, FFTW plans with a progressively wider
// range of algorithms until the the given time limit is reached or the given range of algorithms is explored, returning the best
// available plan.
//
// For example, specifying FFTW_PATIENT first plans in FFTW_ESTIMATE mode, then in FFTW_MEASURE mode, then finally (time permitting)
// in FFTW_PATIENT. If FFTW_EXHAUSTIVE is specified instead, the planner will further progress to FFTW_EXHAUSTIVE mode.
//
// Note that the seconds argument specifies only a rough limit; in practice, the planner may use somewhat more time if the time limit
// is reached when the planner is in the middle of an operation that cannot be interrupted. At the very least, the planner will complete
// planning in FFTW_ESTIMATE mode (which is thus equivalent to a time limit of 0).

const unsigned COOL_DATACACHE_SIZE = 32 * 1048576; // 32 MB

void cool_datacache(unsigned long long memsize)
{
    // We cool the cache by allocating memory and performing the Sieve of Eratosthenes.

    uint8_t * z = new uint8_t[memsize];

    for (unsigned long long i = 0; i < memsize; ++i)
    {
        z[i] = (i & 255);
    }

    delete [] z;
}

unsigned long long cool_datacache_old(unsigned long long memsize)
{
    // We cool the cache by allocating memory and performing the Sieve of Eratosthenes.

    uint8_t * sieve = new uint8_t[memsize];

    for (unsigned long long i = 0; i < memsize; ++i)
    {
        sieve[i] = (i >= 2); // all numbers >= 2 are candidate primes.
    }

    for (unsigned long long i = 2; i < memsize; ++i)
    {
        if (sieve[i] != 0) // promote candidate prime to prime.
        {
            for (unsigned j = i + i; j < memsize; j += i)
            {
                // knock out multiples of i.
                sieve[j] = 0;
            }
        }
    }

    // count primes

    unsigned long long count_primes = 0;

    for (unsigned i = 0; i < memsize; ++i)
    {
        count_primes += sieve[i];
    }

    delete [] sieve;

    return count_primes;
}

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

class Sampler
{
    public:

        Sampler(unsigned n)
        {
            samples.reserve(n);
        }

        void add(const double & x)
        {
            samples.push_back(x);
        }

        double median() const
        {
            const unsigned n = samples.size();

            if (n == 0)
            {
                return std::numeric_limits<double>::quiet_NaN();
            }

            std::vector<double> sorted_samples(samples);
            std::sort(sorted_samples.begin(), sorted_samples.end());

            if (n % 2 != 0)
            {
                // odd number of samples
                return sorted_samples[(n - 1) / 2];
            }
            else
            {
                // even number of samples; return mean of two middle samples.
                return (sorted_samples[n / 2 - 1] + sorted_samples[n / 2]) / 2.0;
            }
        }

        double mean() const
        {
            double sum = accumulate(samples.begin(), samples.end(), 0.0);

            return sum / samples.size();
        }

        double min() const
        {
            return *min_element(samples.begin(), samples.end());
        }

        double max() const
        {
            return *max_element(samples.begin(), samples.end());
        }

        double stddev() const
        {
            const double mu = mean();
            double sum2 = 0.0;

            for (const double & x : samples)
            {
                const double deviation = (x - mu);
                sum2 += deviation * deviation;
            }

            return sqrt(sum2 / (samples.size() - 1));
        }

        unsigned n() const
        {
            return samples.size();
        }

    private:

        std::vector<double> samples;
};

template <typename traits>
void execute_tests_r2c_1d(const std::vector<SignalType> & signal_options, const std::vector<unsigned> & flags_options, const std::vector<CacheTemperature> & cache_options, const unsigned & n_in, const unsigned & repeats)
{
    const unsigned n_out = n_in / 2 + 1;

    for (const SignalType & signalType : signal_options)
    {
        vector<typename traits::real_type> v(n_in);

        switch(signalType)
        {
            case SignalZero: break;
            case SignalDC: fill(v.begin(), v.end(), 1.41592653589793238462); break;
            case SignalNoise:
            {
                unsigned i = 0;
                while (i < v.size())
                {
                    const double U = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
                    const double V = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);

                    const double X = sqrt(-2.0 * log(U)) * cos(2 * M_PI * V);
                    const double Y = sqrt(-2.0 * log(U)) * sin(2 * M_PI * V);

                    if (i < v.size())
                    {
                        v[i++] = X;
                    }

                    if (i < v.size())
                    {
                        v[i++] = Y;
                    }
                }
            }
        }

        for (const unsigned & flags : flags_options)
        {
            typename traits::real_type * x = traits::alloc_real(n_in);
            assert(x != nullptr);

            typename traits::complex_type * y = traits::alloc_complex(n_out);
            assert(y != nullptr);

            auto t1 = chrono::high_resolution_clock::now();
            typename traits::plan plan = traits::plan_dft_r2c_1d(n_in, x, y, flags);
            auto t2 = chrono::high_resolution_clock::now();

            assert(plan != nullptr);

            const double planning_duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count() / 1e6;

            for (const CacheTemperature & cacheTemperature : cache_options)
            {

                // Initialize the input
                // Copy test data

                for (unsigned i = 0; i < n_in; ++i)
                {
                    x[i] = v[i];
                }

                Sampler durationSampler(repeats);

                // Note that we actually do 1 extra evaluation.
                // This allows us to skip the first repeat which may be marred by cache effects.

                if (cacheTemperature == CacheWarm)
                {
                    // We want to do a warm-cache test. So let's heat up the cache by excecuting
                    // the transform three times (a bit overkill, but that doesn't matter).

                    for (unsigned i = 1; i <= 3; ++i)
                    {
                        traits::execute(plan);
                    }
                }

                while (durationSampler.n() < repeats)
                {
                    if (cacheTemperature == CacheCold)
                    {
                        // We want do do a cool-cache test.
                        cool_datacache(COOL_DATACACHE_SIZE);
                    }

                    auto t1 = chrono::high_resolution_clock::now();
                    traits::execute(plan);
                    auto t2 = chrono::high_resolution_clock::now();

                    const double duration = chrono::duration_cast<chrono::microseconds>(t2 - t1).count() / 1e6;

                    durationSampler.add(duration);
                }

                cout << scientific;
                cout.precision(6);

                cout << "operation"         " " << setw( 8) << left  << "r2c"                                        << " "
                        "type"              " " << setw(11) << left  << TypeName<typename traits::real_type>::name() << " "
                        "signal"            " " << setw( 5) << left  << sigtype_to_string(signalType)                << " "
                        "flags"             " " << setw(11) << left  << flags_to_string(flags)                       << " "
                        "cache"             " " << setw( 1) << right << cachetemp_to_string(cacheTemperature)        << " "
                        "n_in"              " " << setw(12) << right << n_in                                         << " "
                        "n_out"             " " << setw(12) << right << n_out                                        << " "
                        "planning"          " " << setw(10) << right << planning_duration                            << " "
                        "reps"              " " << setw( 6) << right << durationSampler.n()                          << " "
                        "min"               " " << setw(12) << right << durationSampler.min()                        << " "
                        "median"            " " << setw(12) << right << durationSampler.median()                     << " "
                        "max"               " " << setw(12) << right << durationSampler.max()                        << " "
                        "mean"              " " << setw(12) << right << durationSampler.mean()                       << " "
                        "stddev"            " " << setw(12) << right << durationSampler.stddev()                     << endl;

                traits::destroy_plan(plan);

                traits::free(y);
                traits::free(x);

                traits::cleanup();

            } // flags loop
        } // cache loop
    } // signal type loop
}

int main()
{
    const unsigned repeats = 10000;

    const std::vector<unsigned> flags_options = {
        //FFTW_ESTIMATE   | FFTW_PRESERVE_INPUT | FFTW_UNALIGNED,
        //FFTW_ESTIMATE   | FFTW_DESTROY_INPUT  | FFTW_UNALIGNED,
        //FFTW_ESTIMATE   | FFTW_PRESERVE_INPUT,
        //FFTW_ESTIMATE   | FFTW_DESTROY_INPUT,
        //FFTW_MEASURE    | FFTW_PRESERVE_INPUT | FFTW_UNALIGNED,
        //FFTW_MEASURE    | FFTW_DESTROY_INPUT  | FFTW_UNALIGNED,
        //FFTW_MEASURE    | FFTW_PRESERVE_INPUT,
        //FFTW_MEASURE    | FFTW_DESTROY_INPUT,
        //FFTW_PATIENT    | FFTW_PRESERVE_INPUT | FFTW_UNALIGNED,
        //FFTW_PATIENT    | FFTW_DESTROY_INPUT  | FFTW_UNALIGNED,
        //FFTW_PATIENT    | FFTW_PRESERVE_INPUT,
        //FFTW_PATIENT    | FFTW_DESTROY_INPUT,
        //FFTW_EXHAUSTIVE | FFTW_PRESERVE_INPUT | FFTW_UNALIGNED,
        //FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT  | FFTW_UNALIGNED,
        //FFTW_EXHAUSTIVE | FFTW_PRESERVE_INPUT,
        FFTW_EXHAUSTIVE | FFTW_DESTROY_INPUT
    };

    const std::vector<CacheTemperature> cache_options = {
        //CacheCold,
        CacheWarm
    };

    const std::vector<SignalType> signal_options = {
        SignalZero,
        SignalDC,
        SignalNoise,
        SignalZero,
        SignalDC,
        SignalNoise
    };

    for (unsigned n = 1; ; n *= 2)
    {
        execute_tests_r2c_1d<FFTW_Traits<double>>(signal_options, flags_options, cache_options, n, repeats);
    }

    return EXIT_SUCCESS;
}
