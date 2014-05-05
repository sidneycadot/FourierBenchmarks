
CXXFLAGS = -W -Wall -O3 -std=c++11 -g
CPPFLAGS = -I..

CC=$(CXX)

#LDLIBS =  # -lfftw3q
LDLIBS = -lfftw3f -lfftw3 -lfftw3l -lmpc -lmpfr -lgmp

#LDLIBS = -lfftw3f_threads -lfftw3f -lfftw3_threads -lfftw3 -lfftw3l_threads -lfftw3l
#LDLIBS = -lfftw3_omp -lfftw3_omp -lfftw3l_omp # -lfftw3q_omp

TARGETS = TestFFTW-performance TestFFTW-precision

vpath % ../ReferenceImplementation
vpath % ../SignalGenerator
vpath % ../MultiPrecisionUtils

.PHONY : clean default

default : $(TARGETS)

TestFFTW-performance : TestFFTW-performance.o FftwUtils.h

TestFFTW-performance.o : TestFFTW-performance.cc FftwUtils.h

TestFFTW-precision : TestFFTW-precision.o ReferenceImplementation.o MultiPrecisionUtils.o

TestFFTW-precision.o : TestFFTW-precision.cc FftwUtils.h SignalGenerator.h ReferenceImplementation.h

ReferenceImplementation.o : ReferenceImplementation.cc ReferenceImplementation.h

MultiPrecisionUtils.o : MultiPrecisionUtils.cc MultiPrecisionUtils.h

TestFFTW : TestFFTW.cc FftwUtils.h

clean :
	$(RM) $(TARGETS) *.o *~
