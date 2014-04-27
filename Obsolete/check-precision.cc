
////////////////////////
// check-precision.cc //
////////////////////////

#include "naive_fourier_transform.h"
#include "recursive_fft.h"
#include "fftw_fft.h"

//#include <cmath>
#include <iostream>
//#include <vector>
//#include <complex>
//#include <cassert>
//#include <chrono>
//#include <fftw3.h>

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

template <typename real_type>
vector<complex<real_type>> test_naive_fourier(const vector<complex<real_type>> & tvec)
{
    vector<complex<real_type>> fvec(tvec.size());

    naive_fourier_transform<real_type>(tvec, fvec);

    return fvec;
}

template <typename real_type>
vector<complex<real_type>> test_recursive_fft(const vector<complex<real_type>> & tvec)
{
    vector<complex<real_type>> fvec = tvec;

    recursive_fft<real_type>(fvec.data(), fvec.size(), 1);

    return fvec;
}

template <typename real_type>
vector<complex<real_type>> test_fftw(const vector<complex<real_type>> & tvec)
{
    vector<complex<real_type>> fvec(tvec.size());

    run_fftw(tvec, fvec);

    return fvec;
}

struct TestResult
{
    TestResult(const string & name, const vector<complex<double>> & result) :
        name(name),
        result(result)
    {
        // empty body
    }

    string name;
    vector<complex<double>> result;
};

void compare(const TestResult & t1, const TestResult & t2)
{
    cout << "=== " << t1.name << " vs. " << t2.name << " ===" << endl;

    const unsigned n = t1.result.size();
    assert(t2.result.size() == n);

    double err = 0;

    for (unsigned i = 0; i < n; ++i)
    {
        double e = abs(t1.result[i] - t2.result[i]);
        err += e * e;
    }

    cout << err << endl;
}

int main()
{
    for (unsigned pow_two = 1; pow_two <= 8; ++pow_two)
    {
        const unsigned n = (1 << pow_two);

        cout << "***** n = " << n << endl << endl;

        vector<complex<double>> tvec = mk_complex_test_vector(n);

        vector<TestResult> results;

        results.push_back(TestResult("naive<double>"     , test_naive_fourier<double>(tvec)));
        results.push_back(TestResult("recursive<double>" , test_recursive_fft<double>(tvec)));
        results.push_back(TestResult("fftw<double>"      , test_fftw<double>(tvec)));

        for (unsigned i = 0; i < results.size(); ++i)
        {
            for (unsigned j = i + 1; j < results.size(); ++j)
            {
                compare(results[i], results[j]);
            }
        }

        cout << endl;
    }

    return 0;
}
