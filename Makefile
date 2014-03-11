
CXXFLAGS = -W -Wall -O3 -std=c++11

LDLIBS = -lfftw3

LDLIBS = -lmkl_core -lmkl_intel_lp64 -lmkl_core -lmkl_intel_lp64 -lmkl_sequential -lpthread

LDFLAGS = -L /opt/intel/composer_xe_2013_sp1.2.144/mkl/lib/intel64

CPPFLAGS = -I /opt/intel/composer_xe_2013_sp1.2.144/mkl/include

check-precision : check-precision.cc

naive : naive.cc

mkl_fft : mkl_fft.cc

clean :
	$(RM) naive check-precision
