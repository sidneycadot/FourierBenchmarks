
///////////////////
// testvectors.h //
///////////////////

#ifndef testvectors_h
#define testvectors_h

#include <vector>
#include <cmath>
#include <complex>

class PrimeGenerator
{
    public:

        static unsigned prime(const unsigned & i)
        {
            while (i >= primevec.size())
            {
                limit *= 2;

                primevec = sieve(limit);
            }
            return primevec[i];
        }

        static std::vector<unsigned> sieve(unsigned int & limit)
        {
            std::vector<unsigned> primevec;

            std::vector<bool> sieve(limit + 1);

            for (unsigned i = 2; i <= limit; ++i)
            {
                if (!sieve[i])
                {
                    primevec.push_back(i);

                    for (unsigned j = i + i; j <= limit; j += i)
                    {
                        sieve[j] = true;
                    }
                }
            }

            return primevec;
        }

    private:

        static unsigned limit;
        static std::vector<unsigned> primevec;
};

template <typename real_type>
std::vector<real_type> mk_real_test_vector(const unsigned & n)
{
    std::vector<real_type> v;

    v.reserve(n);

    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned p = PrimeGenerator::prime(i);

        v.push_back(sqrt(p));
    }

    return v;
}

template <typename real_type>
std::vector<std::complex<real_type>> mk_complex_test_vector(const unsigned & n)
{
    std::vector<std::complex<real_type>> v;

    v.reserve(n);

    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned p = PrimeGenerator::prime(i);

        v.push_back(std::complex<real_type>(sqrt(p), cbrt(p)));
    }

    return v;
}

#endif // testvectors_h
