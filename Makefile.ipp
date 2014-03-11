
.PHONY : clean

IPP_DIR = /opt/intel/composer_xe_2013_sp1.2.144/ipp

CPPFLAGS = -I $(IPP_DIR)/include
LDFLAGS  = -L $(IPP_DIR)/lib/intel64

CXXFLAGS = -std=c++11

LDLIBS = -lippcore -lipps

ipp_fft : ipp_fft.cc


clean :
	$(RM) ipp_fft
