#! /usr/bin/env python

import numpy as np
import time

def nextpow2(n):
    n -= 1
    p = 0
    while n > 0:
        n //= 2
        p += 1
    return p

def matlab_czt(x, k = None, w = None, a = None):

    #CZT  Chirp z-transform.
    #   G = CZT(X, M, W, A) is the M-element z-transform of sequence X,
    #   where M, W and A are scalars which specify the contour in the z-plane
    #   on which the z-transform is computed.  M is the length of the transform,
    #   W is the complex ratio between points on the contour, and A is the
    #   complex starting point.  More explicitly, the contour in the z-plane
    #   (a spiral or "chirp" contour) is described by
    #       z = A * W.^(-(0:M-1))

    #   The parameters M, W, and A are optional; their default values are
    #   M = length(X), W = exp(-1i*2*pi/M), and A = 1.  These defaults
    #   cause CZT to return the z-transform of X at equally spaced points
    #   around the unit circle, equivalent to FFT(X).

    #   If X is a matrix, the chirp z-transform operation is applied to each
    #   column.

    #   See also FFT, FREQZ.

    #   Author(s): C. Denham, 1990.
    #          J. McClellan, 7-25-90, revised
    #          C. Denham, 8-15-90, revised
    #          T. Krauss, 2-16-93, updated help
    #   Copyright 1988-2010 The MathWorks, Inc.
    #       $Revision: 1.7.4.3 $  $Date: 2010/12/06 00:01:36 $

    #   References:
    #     [1] Oppenheim, A.V. & R.W. Schafer, Discrete-Time Signal
    #         Processing,  Prentice-Hall, pp. 623-628, 1989.
    #     [2] Rabiner, L.R. and B. Gold, Theory and Application of
    #         Digital Signal Processing, Prentice-Hall, Englewood
    #         Cliffs, New Jersey, pp. 393-399, 1975.

    m = len(x)

    if k is None:
        k = len(x)

    if w is None:
        w = np.exp(-2j * np.pi / k)

    if a is None:
        a = 1

    # ------- Length for power-of-two fft.

    nfft = 2 ** nextpow2(m + k - 1)

    # ------- Premultiply data.

    kk = np.arange(1 - m, max(k, m))

    ww = w ** ((kk ** 2) / 2.0)                 # <----- Chirp filter is 1/ww

    y = x * a ** -np.arange(m) * ww[m - 1 : 2 * m - 1]

    # ------- Fast convolution via FFT.

    fy = np.fft.fft(y, nfft)

    fv = np.fft.fft(ww[0 : k + m - 1].conj(), nfft)   # <----- Chirp filter.

    fy = fy * fv

    g  = np.fft.ifft(fy)

    # ------- Final multiply.

    g = g[m - 1 : m + k - 1] * ww[m - 1 : m + k - 1]

    return g

def octave_czt(x, m = None, w = None, a = None):

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

    if m is None:
        m = len(x)

    if w is None:
        w = np.exp(-2j * np.pi / m)

    if a is None:
        a = 1

    n = len(x)

    nfft = 2 ** nextpow2(n + m - 1) # fft pad

    W2 = w ** ((np.arange(-(n - 1), max(m, n)) ** 2) / 2.0) # chirp

    fg = np.fft.fft(x * (a ** -np.arange(n)) * W2[n - 1 : 2 * n - 1], nfft)

    fw = np.fft.fft(W2[0 : m + n - 1].conj(), nfft)

    gg = np.fft.ifft(fg * fw, nfft)

    y = gg[n - 1 : m + n - 1] * W2[n - 1 : m + n - 1]

    return y

NUM_POINTS = 999983

x = np.log(123.0 + np.arange(NUM_POINTS)) + np.log(321.0 + np.arange(NUM_POINTS)) * 1j

t1 = time.time()
numpy_fft_result = np.fft.fft(x)
t2 = time.time()
numpy_fft_time = (t2 - t1)
print "numpy_fft_time", numpy_fft_time

t1 = time.time()
matlab_czt_result = matlab_czt(x)
t2 = time.time()
matlab_czt_time = (t2 - t1)

err_matlab_czt = matlab_czt(x) - np.fft.fft(x)
err_matlab_czt = sum(err_matlab_czt.conj() * err_matlab_czt)

print "matlab_czt_time", matlab_czt_time
print "matlab_czt_error", abs(err_matlab_czt)

t1 = time.time()
octave_czt_result = matlab_czt(x)
t2 = time.time()
octave_czt_time = (t2 - t1)

err_octave_czt = octave_czt(x) - np.fft.fft(x)
err_octave_czt = sum(err_octave_czt.conj() * err_octave_czt)

print "octave_czt_time", octave_czt_time
print "octave_czt_error", abs(err_octave_czt)
