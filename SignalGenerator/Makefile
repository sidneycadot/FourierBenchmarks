
CXXFLAGS = -W -Wall -O3 -std=c++11
CC=$(CXX)

.PHONY : clean

LDLIBS = -lmpc -lmpfr -lgmp

vpath % ../MultiPrecisionUtils

TestSignalGenerator : TestSignalGenerator.o MultiPrecisionUtils.o

TestSignalGenerator.o : TestSignalGenerator.cc SignalGenerator.h MultiPrecisionUtils.h

MultiPrecisionUtils.o : MultiPrecisionUtils.cc MultiPrecisionUtils.h

clean :
	$(RM) TestSignalGenerator *.o *~

