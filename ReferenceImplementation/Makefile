
CXXFLAGS = -W -Wall -O3 -std=c++11 -g
CC = $(CXX)
CPPFLAGS = -I..

.PHONY : clean

LDLIBS = -lmpc -lmpfr -lgmp

TARGET = TestReferenceImplementation

OBJS = TestReferenceImplementation.o ReferenceImplementation.o MultiPrecisionUtils.o

$(TARGET) : $(OBJS)

vpath % ../MultiPrecisionUtils

TestReferenceImplementation.o : TestReferenceImplementation.cc ReferenceImplementation.h MultiPrecisionUtils.h

ReferenceImplementation.o : ReferenceImplementation.cc ReferenceImplementation.h MultiPrecisionUtils.h

MultiPrecisionUtils.o : MultiPrecisionUtils.cc MultiPrecisionUtils.h

clean :
	$(RM) $(TARGET) $(OBJS)
