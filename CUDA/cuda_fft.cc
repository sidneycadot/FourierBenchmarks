
#include <cassert>
#include <iostream>
#include <vector>
#include <complex>

#include <cufft.h>
#include <cuda_runtime_api.h>

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
    const unsigned n = 8;

    vector<complex<double>> tvec = mk_complex_test_vector(n);

    cout << "tvec: ";
    print_vector(cout, tvec);
    cout << endl;

    vector<complex<double>> fvec(n);

    cudaError_t cudaErr;

    cufftResult cufft_result;
    cufftHandle plan;

    cufftDoubleComplex * tdata_device;

    cudaErr = cudaMalloc(reinterpret_cast<void **>(&tdata_device), n * sizeof(cufftDoubleComplex));
    assert(cudaErr == cudaSuccess);

    cufftDoubleComplex * fdata_device;

    cudaErr = cudaMalloc(reinterpret_cast<void **>(&fdata_device), n * sizeof(cufftDoubleComplex));
    assert(cudaErr == cudaSuccess);

    cufft_result = cufftPlan1d(&plan, n, CUFFT_Z2Z, 1);
    assert(cufft_result == CUFFT_SUCCESS);

    //

    cudaErr = cudaMemcpy(tdata_device, tvec.data(), n * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice);
    assert(cudaErr == cudaSuccess);

    cufft_result = cufftExecZ2Z(plan, tdata_device, fdata_device, CUFFT_FORWARD);
    assert(cufft_result == CUFFT_SUCCESS);

    cudaErr = cudaDeviceSynchronize();
    assert(cudaErr == cudaSuccess);

    cudaErr = cudaMemcpy(fvec.data(), fdata_device, n * sizeof(cufftDoubleComplex), cudaMemcpyDeviceToHost);
    assert(cudaErr == cudaSuccess);

    cout << "fvec: ";
    print_vector(cout, fvec);
    cout << endl;

    cufft_result = cufftDestroy(plan);
    assert(cufft_result == CUFFT_SUCCESS);

    cudaErr = cudaFree(fdata_device);
    assert(cudaErr == cudaSuccess);

    cudaErr = cudaFree(tdata_device);
    assert(cudaErr == cudaSuccess);

    return 0;
}
