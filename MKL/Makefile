
.PHONY : clean

MKL_DIR = /opt/intel/composer_xe_2013_sp1.2.144/mkl

CXXFLAGS = -std=c++11 -W -Wall -O3
CC=$(CXX)

CPPFLAGS = -I $(MKL_DIR)/include
LDFLAGS  = -L $(MKL_DIR)/lib/intel64

LDLIBS = -lmkl_core -lmkl_intel_ilp64 -lmkl_core -lmkl_sequential -lpthread

TestMKL : TestMKL.o MklDftiUtils.o

TestMKL.o : TestMKL.cc MklDftiUtils.h

MklDftiUtils.o : MklDftiUtils.cc MklDftiUtils.h

clean :
	$(RM) TestMKL TestMKL.o MklDftiUtils.o
