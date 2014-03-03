
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <cassert>

using namespace std;

// Variations:
//
//     Number of transforms
//     Forward vs inverse
//     Scaling
//     Real vs complex
//     Basic data type (float, double, long double)
//
//     Dimensionality of transform
//     Memory layout
//     In-place or not
//
//     Multithreading
//
// Implementations:
//
// CUDA
// MKL
// IPP
// FFTW
// FFTS
// KissFFT
//
// Naive
// Simple
// OpenCL?
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

int main()
{
    vector<complex<double>> v = mk_complex_test_vector(17);

    cout << "v = ";
    print_vector(cout, v);
    cout << ";" << endl;

    vector<complex<double>> fv(v.size());

    slow_fourier_1d(v, fv);

    cout << "fv = ";
    print_vector(cout, fv);
    cout << ";" << endl;

    return 0;
}
