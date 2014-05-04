
////////////////////
// MklDftiUtils.h //
////////////////////

#ifndef MklDftiUtils_h
#define MklDftiUtils_h

#include <string>
#include <ostream>

#include <mkl_dfti.h>

std::string DftiConfigValueToString(const DFTI_CONFIG_VALUE & value);

void print_mkl_dfti_descriptor_info(std::ostream & out, const DFTI_DESCRIPTOR_HANDLE & descriptor);

#endif // MklDftiUtils_h
