
CFLAGS = -W -Wall -O3 -std=c99
LDLIBS = -lm

.PHONY : clean

main : main.o fft_v0.o fft_v1.o czt.o

main.o : main.c fft_v0.h fft_v1.h

fft_v0.o : fft_v0.c fft_v0.h

fft_v1.o : fft_v1.c fft_v1.h

czt.o : czt.c czt.h

clean :
	$(RM) main main.o fft_v0.o fft_v1.o
