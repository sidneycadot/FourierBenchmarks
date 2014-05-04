
/////////////////////
// MklDftiUtils.cc //
/////////////////////

#include <cassert>
#include <ostream>

#include "MklDftiUtils.h"

using namespace std;

string DftiConfigValueToString(const DFTI_CONFIG_VALUE & value)
{
    bool NOT_IMPLEMENTED = false;

    const char * str = nullptr;

    switch (value)
    {
        // DFTI_COMMIT_STATUS

        case DFTI_COMMITTED   : str = "DFTI_COMMITTED"   ; break;
        case DFTI_UNCOMMITTED : str = "DFTI_UNCOMMITTED" ; break;

        // DFTI_FORWARD_DOMAIN

        case DFTI_COMPLEX        : str = "DFTI_COMPLEX"  ; break;
        case DFTI_REAL           : str = "DFTI_REAL"     ; break;

        // case DFTI_CONJUGATE_EVEN : str = "DFTI_CONJUGATE_EVEN"; NOT_IMPLEMENTED = true; break;

        // DFTI_PRECISION

        case DFTI_SINGLE : str = "DFTI_SINGLE"; break;
        case DFTI_DOUBLE : str = "DFTI_DOUBLE"; break;

        // DFTI_FORWARD_SIGN

        // case DFTI_NEGATIVE : str = "DFTI_NEGATIVE"; NOT_IMPLEMENTED = true; break;
        // case DFTI_POSITIVE : str = "DFTI_POSITIVE"; NOT_IMPLEMENTED = true; break;

        // DFTI_COMPLEX_STORAGE and DFTI_CONJUGATE_EVEN_STORAGE

        case DFTI_COMPLEX_COMPLEX : str = "DFTI_COMPLEX_COMPLEX"; break;
        case DFTI_COMPLEX_REAL    : str = "DFTI_COMPLEX_REAL"   ; break;

        // DFTI_REAL_STORAGE
        case DFTI_REAL_COMPLEX : str = "DFTI_REAL_COMPLEX"; break;
        case DFTI_REAL_REAL    : str = "DFTI_REAL_REAL"   ; break;

        // DFTI_PLACEMENT
        case DFTI_INPLACE     : str = "DFTI_INPLACE"    ; break;
        case DFTI_NOT_INPLACE : str = "DFTI_NOT_INPLACE"; break;

        // DFTI_INITIALIZATION_EFFORT

        // case DFTI_LOW    : str = "DFTI_LOW";    NOT_IMPLEMENTED = true; break;
        // case DFTI_MEDIUM : str = "DFTI_MEDIUM"; NOT_IMPLEMENTED = true; break;
        // case DFTI_HIGH   : str = "DFTI_HIGH";   NOT_IMPLEMENTED = true; break;

        // DFTI_ORDERING

        case DFTI_ORDERED            : str = "DFTI_ORDERED"            ; break;
        case DFTI_BACKWARD_SCRAMBLED : str = "DFTI_BACKWARD_SCRAMBLED" ; break;

     // case DFTI_FORWARD_SCRAMBLED  : str = "DFTI_FORWARD_SCRAMBLED"; NOT_IMPLEMENTED = true; break;

        // Allow/avoid certain usages

        case DFTI_ALLOW : str = "DFTI_ALLOW" ; break;
        case DFTI_AVOID : str = "DFTI_AVOID" ; break;
        case DFTI_NONE  : str = "DFTI_NONE"  ; break;

        // DFTI_PACKED_FORMAT (for storing conjugate-even finite sequence in real array)
        case DFTI_CCS_FORMAT  : str = "DFTI_CCS_FORMAT"  ; break;
        case DFTI_PACK_FORMAT : str = "DFTI_PACK_FORMAT" ; break;
        case DFTI_PERM_FORMAT : str = "DFTI_PERM_FORMAT" ; break;
        case DFTI_CCE_FORMAT  : str = "DFTI_CCE_FORMAT"  ; break;
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

void print_mkl_dfti_descriptor_info(ostream & out, const DFTI_DESCRIPTOR_HANDLE & descriptor)
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
}
