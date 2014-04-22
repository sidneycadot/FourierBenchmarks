
CXXFLAGS = -W -Wall -O3 -std=c++11 -g
CC = $(CXX)

.PHONY : clean

LDFLAGS = --static
LDLIBS = -lmpc -lmpfr -lgmp

TARGET = Test_FastFourierTransform_ReferenceImplementation

OBJS = Test_FastFourierTransform_ReferenceImplementation.o FastFourierTransform_ReferenceImplementation.o

$(TARGET) : $(OBJS)

Test_FastFourierTransform_ReferenceImplementation.o : Test_FastFourierTransform_ReferenceImplementation.cc FastFourierTransform_ReferenceImplementation.h

FastFourierTransform_ReferenceImplementation.o : FastFourierTransform_ReferenceImplementation.cc FastFourierTransform_ReferenceImplementation.h

clean :
	$(RM) $(TARGET) $(OBJS)
