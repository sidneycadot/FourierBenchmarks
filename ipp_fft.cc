
#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <ipps.h>
#include <ippcore.h>

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

    //DFTI_DESCRIPTOR_HANDLE descriptor;

    //MKL_LONG status;

    //status = DftiCreateDescriptor(&descriptor, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
    //assert(status == DFTI_NO_ERROR);

    //status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    //assert(status == DFTI_NO_ERROR);

    //status = DftiCommitDescriptor(descriptor);
    //assert(status == DFTI_NO_ERROR);


    Ipp64fc * pSrc = ippsMalloc_64fc(n);
    assert(pSrc != nullptr);

    Ipp64fc * pDst = ippsMalloc_64fc(n);
    assert(pDst != nullptr);

    ippsCopy_64fc(reinterpret_cast<Ipp64fc *>(tvec.data()), pSrc, n);

    // Placeholder for the FFT

    IppStatus status;

    const int flag = IPP_FFT_NODIV_BY_ANY;
    const IppHintAlgorithm hint = ippAlgHintFast;

    int pDftSpecSize, pInitBufferSize, pWorkBufferSize;

    status = ippsDFTGetSize_C_64fc(n, flag, hint, &pDftSpecSize, &pInitBufferSize, &pWorkBufferSize);
    assert(status == ippStsNoErr);

    cout << "pDftSpecSize: "    << pDftSpecSize    << endl;
    cout << "pInitBufferSize: " << pInitBufferSize << endl;
    cout << "pWorkBufferSize: " << pWorkBufferSize << endl;

    IppsDFTSpec_C_64fc * pDftSpec = reinterpret_cast<IppsDFTSpec_C_64fc *>(ippsMalloc_8u(pDftSpecSize));
    assert(pDftSpec != nullptr);

    Ipp8u * pInitBuffer = nullptr;
    if (pInitBufferSize != 0)
    {
        pInitBuffer = ippsMalloc_8u(pInitBufferSize);
        assert(pInitBuffer != nullptr);
    }

    Ipp8u * pWorkBuffer = nullptr;
    if (pWorkBufferSize != 0)
    {
        pWorkBuffer = ippsMalloc_8u(pWorkBufferSize);
        assert(pWorkBuffer != nullptr);
    }

    status = ippsDFTInit_C_64fc(n, flag, hint, pDftSpec, pInitBuffer);
    cout << "status: " << ippGetStatusString(status) << endl;
    assert(status == ippStsNoErr);

    if (pInitBuffer != nullptr)
    {
        ippsFree(pInitBuffer);
    }

    status = ippsDFTFwd_CToC_64fc(pSrc, pDst, pDftSpec, pWorkBuffer);
    cout << "status: " << ippGetStatusString(status) << endl;
    assert(status == ippStsNoErr);

    vector<complex<double>> fvec(n);

    ippsCopy_64fc(pDst, reinterpret_cast<Ipp64fc *>(fvec.data()), n);

    cout << "fvec: ";
    print_vector(cout, fvec);
    cout << endl;

    ippsFree(pDftSpec);

    if (pWorkBuffer != nullptr)
    {
        ippsFree(pWorkBuffer);
    }

    ippsFree(pDst);
    ippsFree(pSrc);

    return 0;
}
