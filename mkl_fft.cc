
#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <mkl_dfti.h>
#include <mkl_service.h>

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

string DftiConfigValueToString(const DFTI_CONFIG_VALUE & value)
{
    bool NOT_IMPLEMENTED = false;
    const char * str = nullptr;
    switch (value)
    {
        // DFTI_COMMIT_STATUS
        case DFTI_COMMITTED   : str = "DFTI_COMMITTED"; break;
        case DFTI_UNCOMMITTED : str = "DFTI_UNCOMMITTED"; break;

        // DFTI_FORWARD_DOMAIN
        case DFTI_COMPLEX        : str = "DFTI_COMPLEX"; break;
        case DFTI_REAL           : str = "DFTI_REAL"; break;
     // case DFTI_CONJUGATE_EVEN : str = "DFTI_CONJUGATE_EVEN"; NOT_IMPLEMENTED = true; break;

        // DFTI_PRECISION
        case DFTI_SINGLE : str = "DFTI_SINGLE"; break;
        case DFTI_DOUBLE : str = "DFTI_DOUBLE"; break;

        // DFTI_FORWARD_SIGN
     // case DFTI_NEGATIVE : str = "DFTI_NEGATIVE"; NOT_IMPLEMENTED = true; break;
     // case DFTI_POSITIVE : str = "DFTI_POSITIVE"; NOT_IMPLEMENTED = true; break;

        // DFTI_COMPLEX_STORAGE and DFTI_CONJUGATE_EVEN_STORAGE

        case DFTI_COMPLEX_COMPLEX : str = "DFTI_COMPLEX_COMPLEX"; break;
        case DFTI_COMPLEX_REAL    : str = "DFTI_COMPLEX_REAL"; break;

        // DFTI_REAL_STORAGE
        case DFTI_REAL_COMPLEX : str = "DFTI_REAL_COMPLEX"; break;
        case DFTI_REAL_REAL    : str = "DFTI_REAL_REAL"; break;

        // DFTI_PLACEMENT
        case DFTI_INPLACE     : str = "DFTI_INPLACE"; break;
        case DFTI_NOT_INPLACE : str = "DFTI_NOT_INPLACE"; break;

        // DFTI_INITIALIZATION_EFFORT
     // case DFTI_LOW    : str = "DFTI_LOW";    NOT_IMPLEMENTED = true; break;
     // case DFTI_MEDIUM : str = "DFTI_MEDIUM"; NOT_IMPLEMENTED = true; break;
     // case DFTI_HIGH   : str = "DFTI_HIGH";   NOT_IMPLEMENTED = true; break;

        // DFTI_ORDERING
        case DFTI_ORDERED            : str = "DFTI_ORDERED"; break;
        case DFTI_BACKWARD_SCRAMBLED : str = "DFTI_BACKWARD_SCRAMBLED"; break;
     // case DFTI_FORWARD_SCRAMBLED  : str = "DFTI_FORWARD_SCRAMBLED"; NOT_IMPLEMENTED = true; break;

        // Allow/avoid certain usages
        case DFTI_ALLOW : str = "DFTI_ALLOW"; break;
        case DFTI_AVOID : str = "DFTI_AVOID"; break;
        case DFTI_NONE  : str = "DFTI_NONE"; break;

        // DFTI_PACKED_FORMAT (for storing conjugate-even finite sequence in real array)
        case DFTI_CCS_FORMAT  : str = "DFTI_CCS_FORMAT"; break;
        case DFTI_PACK_FORMAT : str = "DFTI_PACK_FORMAT"; break;
        case DFTI_PERM_FORMAT : str = "DFTI_PERM_FORMAT"; break;
        case DFTI_CCE_FORMAT  : str = "DFTI_CCE_FORMAT"; break;
    }

    if (str)
    {
        if (NOT_IMPLEMENTED)
        {
            return string(str) + " (" + to_string(value) + ") --- NOT IMPLEMENTED";
        }
        else
        {
            return string(str) + " (" + to_string(value) + ")";
        }
    }
    else
    {
        return "unknown value (" + to_string(value) + ")";
    }
}

#if 0
enum DFTI_CONFIG_PARAM
{
    /* Domain for forward transform. No default value */
    DFTI_FORWARD_DOMAIN = 0,

    /* Dimensionality, or rank. No default value */
    DFTI_DIMENSION = 1,

    /* Length(s) of transform. No default value */
    DFTI_LENGTHS = 2,

    /* Floating point precision. No default value */
    DFTI_PRECISION = 3,

    /* Scale factor for forward transform [1.0] */
    DFTI_FORWARD_SCALE  = 4,

    /* Scale factor for backward transform [1.0] */
    DFTI_BACKWARD_SCALE = 5,

    /* Exponent sign for forward transform [DFTI_NEGATIVE]  */
    /* DFTI_FORWARD_SIGN = 6, ## NOT IMPLEMENTED */

    /* Number of data sets to be transformed [1] */
    DFTI_NUMBER_OF_TRANSFORMS = 7,

    /* Storage of finite complex-valued sequences in complex domain
       [DFTI_COMPLEX_COMPLEX] */
    DFTI_COMPLEX_STORAGE = 8,

    /* Storage of finite real-valued sequences in real domain
       [DFTI_REAL_REAL] */
    DFTI_REAL_STORAGE = 9,

    /* Storage of finite complex-valued sequences in conjugate-even
       domain [DFTI_COMPLEX_REAL] */
    DFTI_CONJUGATE_EVEN_STORAGE = 10,

    /* Placement of result [DFTI_INPLACE] */
    DFTI_PLACEMENT = 11,

    /* Generalized strides for input data layout [tigth, row-major for
       C] */
    DFTI_INPUT_STRIDES = 12,

    /* Generalized strides for output data layout [tight, row-major
       for C] */
    DFTI_OUTPUT_STRIDES = 13,

    /* Distance between first input elements for multiple transforms
       [0] */
    DFTI_INPUT_DISTANCE = 14,

    /* Distance between first output elements for multiple transforms
       [0] */
    DFTI_OUTPUT_DISTANCE = 15,

    /* Effort spent in initialization [DFTI_MEDIUM] */
    /* DFTI_INITIALIZATION_EFFORT = 16, ## NOT IMPLEMENTED */

    /* Use of workspace during computation [DFTI_ALLOW] */
    DFTI_WORKSPACE = 17,

    /* Ordering of the result [DFTI_ORDERED] */
    DFTI_ORDERING = 18,

    /* Possible transposition of result [DFTI_NONE] */
    DFTI_TRANSPOSE = 19,

    /* User-settable descriptor name [""] */
    DFTI_DESCRIPTOR_NAME = 20, /* DEPRECATED */

    /* Packing format for DFTI_COMPLEX_REAL storage of finite
       conjugate-even sequences [DFTI_CCS_FORMAT] */
    DFTI_PACKED_FORMAT = 21,

    /* Commit status of the descriptor - R/O parameter */
    DFTI_COMMIT_STATUS = 22,

    /* Version string for this DFTI implementation - R/O parameter */
    DFTI_VERSION = 23,

    /* Ordering of the forward transform - R/O parameter */
    /* DFTI_FORWARD_ORDERING  = 24, ## NOT IMPLEMENTED */

    /* Ordering of the backward transform - R/O parameter */
    /* DFTI_BACKWARD_ORDERING = 25, ## NOT IMPLEMENTED */

    /* Number of user threads that share the descriptor [1] */
    DFTI_NUMBER_OF_USER_THREADS = 26,

    /* Limit the number of threads used by this descriptor [0 = don't care] */
    DFTI_THREAD_LIMIT = 27
};
#endif

ostream & show_mkl_dfti_descriptor_info(ostream & out, const DFTI_DESCRIPTOR_HANDLE & descriptor)
{
    // DFTI_PRECISION

    {
        DFTI_CONFIG_VALUE precision;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_PRECISION, &precision);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_PRECISION ................... : " << DftiConfigValueToString(precision) << endl;
    }

    // DFTI_FORWARD_DOMAIN

    {
        DFTI_CONFIG_VALUE forward_domain;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_FORWARD_DOMAIN, &forward_domain);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_FORWARD_DOMAIN .............. : " << DftiConfigValueToString(forward_domain) << endl;
    }

    // DFTI_DIMENSION

    MKL_LONG number_of_dimensions; // number of dimensions

    {
        MKL_LONG status = DftiGetValue(descriptor, DFTI_DIMENSION, &number_of_dimensions);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_DIMENSION ................... : " << number_of_dimensions << endl;
    }

    // DFTI_LENGTH

    {
        MKL_LONG dim[number_of_dimensions];

        MKL_LONG status = DftiGetValue(descriptor, DFTI_LENGTHS, dim);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_LENGTHS ..................... : {";

        for (unsigned i = 0; i < number_of_dimensions; ++i)
        {
            if (i != 0)
            {
                out << ", ";
            }
            out << dim[i];
        }

        out << "}" << endl;
    }

    // DFTI_PLACEMENT

    {
        DFTI_CONFIG_VALUE placement;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_PLACEMENT, &placement);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_PLACEMENT ................... : " << DftiConfigValueToString(placement) << endl;
    }

    // DFTI_FORWARD_SCALE

    {
        double forward_scale;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_FORWARD_SCALE, &forward_scale);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_FORWARD_SCALE ............... : " << forward_scale << endl;
    }

    // DFTI_BACKWARD_SCALE

    {
        double backward_scale;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_FORWARD_SCALE, &backward_scale);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_BACKWARD_SCALE .............. : " << backward_scale << endl;
    }

    // DFTI_NUMBER_OF_USER_THREADS

    {
        MKL_LONG number_of_user_threads;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_NUMBER_OF_USER_THREADS, &number_of_user_threads);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_NUMBER_OF_USER_THREADS ...... : " << number_of_user_threads << endl;
    }

    // DFTI_THREAD_LIMIT

    {
        MKL_LONG thread_limit;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_THREAD_LIMIT, &thread_limit);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_THREAD_LIMIT ................ : " << thread_limit << endl;
    }

    // DFTI_DESCRIPTOR_NAME

    {
        // Actual max length is (DFTI_MAX_NAME_LENGTH - 1); also reserve room for terminating NUL character.
        char descriptor_name[DFTI_MAX_NAME_LENGTH];

        MKL_LONG status = DftiGetValue(descriptor, DFTI_DESCRIPTOR_NAME, &descriptor_name);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_DESCRIPTOR_NAME ............. : \"" << descriptor_name << "\"" << endl;
    }

    // DFTI_INPUT_STRIDES

    {
        MKL_LONG input_strides[number_of_dimensions];

        MKL_LONG status = DftiGetValue(descriptor, DFTI_INPUT_STRIDES, input_strides);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_INPUT_STRIDES ............... : {";

        for (unsigned i = 0; i < number_of_dimensions; ++i)
        {
            if (i != 0)
            {
                out << ", ";
            }
            out << input_strides[i];
        }

        out << "}" << endl;
    }

    // DFTI_OUTPUT_STRIDES

    {
        MKL_LONG output_strides[number_of_dimensions];

        MKL_LONG status = DftiGetValue(descriptor, DFTI_OUTPUT_STRIDES, output_strides);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_OUTPUT_STRIDES .............. : {";

        for (unsigned i = 0; i < number_of_dimensions; ++i)
        {
            if (i != 0)
            {
                out << ", ";
            }
            out << output_strides[i];
        }

        out << "}" << endl;
    }

    // DFTI_NUMBER_OF_TRANSFORMS

    {
        MKL_LONG number_of_transforms;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_NUMBER_OF_TRANSFORMS, &number_of_transforms);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_NUMBER_OF_TRANSFORMS ........ : " << number_of_transforms << endl;
    }

    // DFTI_INPUT_DISTANCE

    {
        MKL_LONG input_distance;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_INPUT_DISTANCE, &input_distance);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_INPUT_DISTANCE .............. : " << input_distance << endl;
    }

    // DFTI_OUTPUT_DISTANCE

    {
        MKL_LONG output_distance;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_OUTPUT_DISTANCE, &output_distance);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_OUTPUT_DISTANCE ............. : " << output_distance << endl;
    }

    // DFTI_COMPLEX_STORAGE

    {
        DFTI_CONFIG_VALUE complex_storage;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_COMPLEX_STORAGE, &complex_storage);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_COMPLEX_STORAGE ............. : " << DftiConfigValueToString(complex_storage) << endl;
    }

    // DFTI_REAL_STORAGE

    {
        DFTI_CONFIG_VALUE real_storage;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_REAL_STORAGE, &real_storage);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_REAL_STORAGE ................ : " << DftiConfigValueToString(real_storage) << endl;
    }

    // DFTI_CONJUGATE_EVEN_STORAGE

    {
        DFTI_CONFIG_VALUE conjugate_even_storage;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_CONJUGATE_EVEN_STORAGE, &conjugate_even_storage);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_CONJUGATE_EVEN_STORAGE ...... : " << DftiConfigValueToString(conjugate_even_storage) << endl;
    }

    // DFTI_PACKED_FORMAT

    {
        DFTI_CONFIG_VALUE packed_format;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_PACKED_FORMAT, &packed_format);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_PACKED_FORMAT ............... : " << DftiConfigValueToString(packed_format) << endl;
    }

    // DFTI_WORKSPACE

    {
        DFTI_CONFIG_VALUE workspace;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_WORKSPACE, &workspace);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_WORKSPACE ................... : " << DftiConfigValueToString(workspace) << endl;
    }

    // DFTI_ORDERING

    {
        DFTI_CONFIG_VALUE ordering;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_ORDERING, &ordering);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_ORDERING .................... : " << DftiConfigValueToString(ordering) << endl;
    }

    // DFTI_COMMIT_STATUS

    {
        DFTI_CONFIG_VALUE commit_status;

        MKL_LONG status = DftiGetValue(descriptor, DFTI_COMMIT_STATUS, &commit_status);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_COMMIT_STATUS ............... : " << DftiConfigValueToString(commit_status) << endl;
    }

    // DFTI_VERSION

    {
        // Actual max length is (DFTI_VERSION_LENGTH - 1); also reserve room for terminating NUL character.
        char dft_version[DFTI_VERSION_LENGTH];

        MKL_LONG status = DftiGetValue(descriptor, DFTI_VERSION, &dft_version);
        assert(status == DFTI_NO_ERROR);

        out << "DFTI_VERSION ..................... : \"" << dft_version << "\"" << endl;
    }

    return out;
}

int main()
{
    const unsigned n = 8;

    vector<complex<double>> tvec = mk_complex_test_vector(n);

    cout << "tvec: ";
    print_vector(cout, tvec);
    cout << endl;

    DFTI_DESCRIPTOR_HANDLE descriptor;

    MKL_LONG status;

    status = DftiCreateDescriptor(&descriptor, DFTI_DOUBLE, DFTI_COMPLEX, 1, n);
    assert(status == DFTI_NO_ERROR);

    status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    assert(status == DFTI_NO_ERROR);

    status = DftiCommitDescriptor(descriptor);
    assert(status == DFTI_NO_ERROR);

    show_mkl_dfti_descriptor_info(cout, descriptor);

    vector<complex<double>> fvec(n);

    status = DftiComputeForward(descriptor, tvec.data(), fvec.data());
    assert(status == DFTI_NO_ERROR);

    status = DftiFreeDescriptor(&descriptor);
    assert(status == DFTI_NO_ERROR);

    cout << "fvec: ";
    print_vector(cout, fvec);
    cout << endl;

    mkl_free_buffers();

    return 0;
}
