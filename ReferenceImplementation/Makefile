
CXXFLAGS = -W -Wall -O3 -std=c++11 -g
CC = $(CXX)
CPPFLAGS = -I../MultiPrecisionUtils

.PHONY : clean

LDFLAGS = --static
LDLIBS = -lmpc -lmpfr -lgmp

TARGET = Test_FastFourierTransform_ReferenceImplementation

OBJS = Test_FastFourierTransform_ReferenceImplementation.o FastFourierTransform_ReferenceImplementation.o MultiPrecisionUtils.o

$(TARGET) : $(OBJS)

vpath % ../MultiPrecisionUtils

Test_FastFourierTransform_ReferenceImplementation.o : Test_FastFourierTransform_ReferenceImplementation.cc FastFourierTransform_ReferenceImplementation.h MultiPrecisionUtils.h

FastFourierTransform_ReferenceImplementation.o : FastFourierTransform_ReferenceImplementation.cc FastFourierTransform_ReferenceImplementation.h MultiPrecisionUtils.h

MultiPrecisionUtils.o : MultiPrecisionUtils.cc MultiPrecisionUtils.h

clean :
	$(RM) $(TARGET) $(OBJS)
