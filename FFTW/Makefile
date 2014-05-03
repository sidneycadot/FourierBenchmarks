
CXXFLAGS = -W -Wall -O3 -std=c++11 -g

#LDLIBS =  # -lfftw3q
LDLIBS = -lfftw3f -lfftw3 -lfftw3l
#LDLIBS = -lfftw3f_threads -lfftw3f -lfftw3_threads -lfftw3 -lfftw3l_threads -lfftw3l
#LDLIBS = -lfftw3_omp -lfftw3_omp -lfftw3l_omp # -lfftw3q_omp

.PHONY : clean

TestFFTW : TestFFTW.cc

clean :
	$(RM) TestFFTW
