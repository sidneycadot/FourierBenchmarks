
CXXFLAGS = -W -Wall -O3 -std=c++11 -g
CC = $(CXX)

.PHONY : clean

LDLIBS = -lmpc

TestReferenceImplementation : TestReferenceImplementation.o ReferenceImplementation.o

TestReferenceImplementation.o : TestReferenceImplementation.cc ReferenceImplementation.h
ReferenceImplementation.o     : ReferenceImplementation.cc ReferenceImplementation.h

clean :
	$(RM) TestReferenceImplementation
