
CXXFLAGS = -W -Wall -O3 -std=c++11

LDLIBS = -lfftw3

naive : naive.cc

clean :
	$(RM) naive
