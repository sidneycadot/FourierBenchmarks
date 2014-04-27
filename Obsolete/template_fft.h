
////////////////////
// template_fft.h //
////////////////////

#ifndef template_fft_h
#define template_fft_h

template <typename real_type, unsigned n, unsigned stride>
struct template_fft_struct
{
    static void run(complex<real_type> * z)
    {
        static_assert(n > 1, "oops");

        template_fft_struct<n / 2, 2 * stride>::run(z         );
        template_fft_struct<n / 2, 2 * stride>::run(z + stride);

        // Now, we need to combine the values of the sub-FFTs

        complex<real_type> temp[n];

        for (unsigned i = 0; i < n / 2; ++i)
        {
            const real_type turn = -2.0 * M_PI * i / n;

            const complex<real_type> coeff = polar(1.0, turn);

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

template <typename real_type, unsigned stride>
struct template_fft_struct<1, stride>
{
    static void run(complex<real_type> * z)
    {
        (void)z; // no-op
    }
};

template <typename real_type, unsigned n>
void template_fft(complex<real_type> * z)
{
    template_fft_struct<n, 1>::run(z);
}

#endif // template_fft_h
