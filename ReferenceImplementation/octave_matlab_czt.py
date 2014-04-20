#! /usr/bin/env python

import numpy as np
import time

def czt(x, m = None, w = None, a = None):

    # Copyright (C) 2004 Daniel Gunyan
    #
    # This program is free software; you can redistribute it and/or modify it under
    # the terms of the GNU General Public License as published by the Free Software
    # Foundation; either version 3 of the License, or (at your option) any later
    # version.
    #
    # This program is distributed in the hope that it will be useful, but WITHOUT
    # ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    # FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
    # details.
    #
    # You should have received a copy of the GNU General Public License along with
    # this program; if not, see <http://www.gnu.org/licenses/>.

    # -*- texinfo -*-
    # @deftypefn  {Function File} {} czt (@var{x})
    # @deftypefnx {Function File} {} czt (@var{x}, @var{m})
    # @deftypefnx {Function File} {} czt (@var{x}, @var{m}, @var{w})
    # @deftypefnx {Function File} {} czt (@var{x}, @var{m}, @var{w}, @var{a})
    # Chirp z-transform.  Compute the frequency response starting at a and
    # stepping by w for m steps.  a is a point in the complex plane, and
    # w is the ratio between points in each step (i.e., radius increases
    # exponentially, and angle increases linearly).
    #
    # To evaluate the frequency response for the range f1 to f2 in a signal
    # with sampling frequency Fs, use the following:
    #
    # @example
    # @group
    # m = 32;                               ## number of points desired
    # w = exp(-j*2*pi*(f2-f1)/((m-1)*Fs));  ## freq. step of f2-f1/m
    # a = exp(j*2*pi*f1/Fs);                ## starting at frequency f1
    # y = czt(x, m, w, a);
    # @end group
    # @end example
    #
    # If you don't specify them, then the parameters default to a fourier
    # transform:
    #     m=length(x), w=exp(-j*2*pi/m), a=1
    #
    # If x is a matrix, the transform will be performed column-by-column.
    # @end deftypefn

    # Algorithm (based on Oppenheim and Schafer, "Discrete-Time Signal
    # Processing", pp. 623-628):
    #   make chirp of length -N+1 to max(N-1,M-1)
    #     chirp => w^([-N+1:max(N-1,M-1)]^2/2)
    #   multiply x by chirped a and by N-elements of chirp, and call it g
    #   convolve g with inverse chirp, and call it gg
    #     pad ffts so that multiplication works
    #     ifft(fft(g)*fft(1/chirp))
    #   multiply gg by M-elements of chirp and call it done

    n = len(x)

    if m is None: m = n
    if w is None: w = np.exp(-2j * np.pi / m)
    if a is None: a = 1

    wExponents1 =   np.arange(    0   , n) ** 2 / 2.0    # n elements         [ 0 .. (n - 1) ]
    wExponents2 = - np.arange(-(n - 1), m) ** 2 / 2.0    # m + n - 1 elements [ -(n - 1) .. +(m - 1) ]
    wExponents3 =   np.arange(    0   , m) ** 2 / 2.0    # m elements         [ 0        ..  (m - 1) ]

    xx = x * a ** -np.arange(n) * w ** wExponents1

    # Determine next-biggest FFT of power-of-two

    nfft = 1
    while nfft < (m + n - 1):
        nfft += nfft

    fxx = np.fft.fft(xx, nfft)

    ww = w ** wExponents2

    fww = np.fft.fft(ww, nfft)

    fyy = fxx * fww

    yy = np.fft.ifft(fyy, nfft)

    # select output

    yy = yy[n - 1 : m + n - 1]

    y = yy * w ** wExponents3

    return y

NUM_POINTS = 3

x = (10.0 + np.arange(NUM_POINTS)) + (20.0 + np.arange(NUM_POINTS) ** 2) * 1j

t1 = time.time()
numpy_fft_result = np.fft.fft(x)
t2 = time.time()
numpy_fft_time = (t2 - t1)
print "numpy_fft_time", numpy_fft_time

t1 = time.time()
czt_result = czt(x)
t2 = time.time()
czt_time = (t2 - t1)

err_czt = czt_result - np.fft.fft(x)
err_czt = sum(err_czt.conj() * err_czt)

print "czt_time", czt_time
print "czt_error", abs(err_czt)
