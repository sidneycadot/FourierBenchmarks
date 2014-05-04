
CXX = nvcc
LDLIBS = -lcufft
CXXFLAGS = --gpu-architecture=sm_13
.PHIONY : clean

TestCUDA : TestCUDA.cc

clean :
	$(RM) TestCUDA
