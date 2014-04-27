
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <chrono>
#include <fftw3.h>

using namespace std;

bool is_prime(unsigned n)
{
    if (n < 2)
    {
        return false;
    }

    for (unsigned d = 2; d * d <= n; ++d)
    {
        if (n % d == 0)
        {
            return false;
        }
    }

    return true;
}

vector<double> mk_real_test_vector(unsigned n)
{
    vector<double> v;
    v.reserve(n);

    unsigned p = 0;
    while (v.size() < n)
    {
        while (!is_prime(p))
        {
            ++p;
        }

        v.push_back(sqrt(p));

        ++p;
    }

    return v;
}

vector<complex<double>> mk_complex_test_vector(unsigned n)
{
    vector<complex<double>> v;
    v.reserve(n);

    unsigned p = 0;
    while (v.size() < n)
    {
        while (!is_prime(p))
        {
            ++p;
        }

        v.push_back(complex<double>(sqrt(p), cbrt(p)));

        ++p;
    }

    return v;
}

void slow_fourier_1d(const vector<complex<double>> & v, vector<complex<double>> & fv)
{
    const unsigned n = v.size();

    // Make sure the destination vector has the correct size.
    assert(fv.size() == n);

    for (unsigned fi = 0; fi < n; ++fi)
    {
        fv[fi] = 0;

        for (unsigned ti = 0; ti < n; ++ti)
        {
            fv[fi] += polar(1.0, -2.0 * M_PI * fi * ti / n) * v[ti];
        }
    }
}

void slow_inverse_fourier_1d(const vector<complex<double>> & fv, vector<complex<double>> & ifv)
{
    const unsigned n = fv.size();

    // Make sure the destination vector has the correct size.
    assert(ifv.size() == n);

    for (unsigned ti = 0; ti < n; ++ti)
    {
        ifv[ti] = 0;

        for (unsigned fi = 0; fi < n; ++fi)
        {
            ifv[ti] += polar(1.0, +2.0 * M_PI * fi * ti / n) * fv[fi] / static_cast<double>(n);
        }
    }
}

void fft(complex<double> * z, unsigned n, unsigned stride)
{
    if (n == 1)
    {
        return;
    }

    assert (n % 2 == 0);

    fft(z         , n / 2, 2 * stride);
    fft(z + stride, n / 2, 2 * stride);

    // Now, we need to combine the values of the sub-FFTs

    complex<double> temp[n];

    for (unsigned i = 0; i < n / 2; ++i)
    {
        const double turn = -2.0 * M_PI * i / n;

        const complex<double> coeff = polar(1.0, turn);

        temp[2 * i + 0] = z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)];
        temp[2 * i + 1] = z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)];
    }

    // And re-arrange them ...

    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned j = (i / 2) + (i % 2) * (n / 2);

        assert(j < n);

        z[j * stride] = temp[i];
    }
}

template <unsigned n, unsigned stride>
struct fft_template_struct
{
    static void run(complex<double> * z)
    {
        static_assert(n > 1, "oops");

        fft_template_struct<n / 2, 2 * stride>::run(z         );
        fft_template_struct<n / 2, 2 * stride>::run(z + stride);

        // Now, we need to combine the values of the sub-FFTs

        complex<double> temp[n];

        for (unsigned i = 0; i < n / 2; ++i)
        {
            const double turn = -2.0 * M_PI * i / n;

            const complex<double> coeff = polar(1.0, turn);

            temp[2 * i + 0] = z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)];
            temp[2 * i + 1] = z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)];
        }

        // And re-arrange them ...

        for (unsigned i = 0; i < n; ++i)
        {
            const unsigned j = (i / 2) + (i % 2) * (n / 2);

            assert(j < n);

            z[j * stride] = temp[i];
        }
    }
};

template <unsigned stride>
struct fft_template_struct<1, stride>
{
    static void run(complex<double> * z)
    {
        (void)z; // no-op
    }
};

template <unsigned n>
void fft(complex<double> * z)
{
    fft_template_struct<n, 1>::run(z);
}

template <typename T>
ostream & print_vector(ostream & os, const vector<T> & v)
{
    os << "{";
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (i != 0)
        {
            os << ", ";
        }
        os << v[i];
    }
    os << "}";

    return os;
}

double test_naive(const unsigned n, const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    vector<complex<double>> fv(n);

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        slow_fourier_1d(v, fv);
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    return duration / repeats;
}

double test_fftw(const unsigned n, const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    fftw_complex * x = fftw_alloc_complex(n);
    assert(x != nullptr);

    fftw_complex * y = fftw_alloc_complex(n);
    assert(y != nullptr);

    fftw_plan plan = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    assert(plan != nullptr);

    // Initialize the input

    // copy test data
    for (unsigned i = 0; i < n; ++i)
    {
        x[i][0] = real(v[i]);
        x[i][1] = imag(v[i]);
    }

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        fftw_execute(plan);
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    fftw_destroy_plan(plan);

    fftw_free(y);
    fftw_free(x);

    fftw_cleanup();

    return duration / repeats;
}

double test_simple(const unsigned n, const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    vector<complex<double>> fv = mk_complex_test_vector(n);

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        fv = v;
        fft(fv.data(), fv.size(), 1);
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    return duration / repeats;
}

template <unsigned n>
double test_template(const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    vector<complex<double>> fv = mk_complex_test_vector(n);

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        fv = v;
        fft<n>(fv.data());
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    return duration / repeats;
}

void show1()
{
    const unsigned n = 8;

    vector<complex<double>> v = mk_complex_test_vector(n);

    print_vector(cout, v); cout << endl;

    vector<complex<double>> fv = v;

    slow_inverse_fourier_1d(v, fv);

    print_vector(cout, fv); cout << endl;
}

void show2()
{
    const unsigned n = 8;

    vector<complex<double>> v = mk_complex_test_vector(n);

    print_vector(cout, v); cout << endl;

    vector<complex<double>> fv = v;

    //slow_fourier_1d(v, fv);
    fft<n>(fv.data());

    print_vector(cout, fv); cout << endl;
}

int main()
{
    //show1();
    //show2();
    //return 0;

    if (0)
    for (unsigned n = 1; n <= (1u<<20); n *= 2)
    {
        double duration_fftw   = test_fftw(n, 10);
        double duration_simple = test_simple(n, 10);

        double speedup = duration_simple / duration_fftw;

        cout << n << "\t" << duration_fftw << "\t" << duration_simple << "\t" << speedup << endl;
    }

    const unsigned n = 65536;

    double duration_fftw   = test_fftw(n, 100);
    double duration_simple = test_template<n>(100);

    double speedup = duration_simple / duration_fftw;

    cout << n << "\t" << duration_fftw << "\t" << duration_simple << "\t" << speedup << endl;

    return 0;
}
