#!/usr/bin/env python

#*********************************************************************************
# DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
#
# This material is based upon work supported under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the U.S. Air Force.
#
# (c) 2023 Massachusetts Institute of Technology.
#
# Subject to FAR52.227-11 Patent Rights - Ownership by the contractor (May 2014)
#
# The software/firmware is provided to you on an As-Is basis
#
# Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
#********************************************************************************/

# Benchmarks cross-correlation to see which is the fastest
# (native, via convolution, or via FFT)

import numpy as np
import time
import pyfftw

def xcorr(a, b):

    # Works like xcorr() in MATLAB

    # Zero pad one of the inputs so they are the same length
    n_pad = len(a) - len(b)
    if n_pad > 0:
        b = np.concatenate([b, np.zeros(n_pad)])
    elif n_pad < 0:
        a = np.concatenate([a, np.zeros(-n_pad)])

    # Perform the cross correlation
    r = np.correlate(a, b, 'full')

    # Lags
    max_lag = len(a) - 1
    l = np.arange(-max_lag, max_lag+1)

    # Done
    return r, l

def xcorr_via_conv(a, b):

    # Works like xcorr() in MATLAB
    # Uses convolution instead of correlation

    # Zero pad one of the inputs so they are the same length
    n_pad = len(a) - len(b)
    if n_pad > 0:
        b = np.concatenate([b, np.zeros(n_pad)])
    elif n_pad < 0:
        a = np.concatenate([a, np.zeros(-n_pad)])

    # Perform the cross correlation via convolution
    r = np.convolve(a, np.conj(b[::-1]))

    # Lags
    max_lag = len(a) - 1
    l = np.arange(-max_lag, max_lag+1)

    # Done
    return r, l

def xcorr_via_fft(a, b):

    # Works like xcorr() in MATLAB
    # Uses FFT for convolution

    # Zero pad one of the inputs so they are the same length
    n_pad = len(a) - len(b)
    if n_pad > 0:
        b = np.concatenate([b, np.zeros(n_pad)])
    elif n_pad < 0:
        a = np.concatenate([a, np.zeros(-n_pad)])

    # Perform forward FFT, multiply, and inverse FFT
    # Note that we need to zero pad since multiplication of the FFTs is circular
    # convolution but we want linear convolution
    n_pad = len(a) - 1
    a_f = np.fft.fft(np.concatenate([a, np.zeros(n_pad)]))
    b_f = np.fft.fft(np.concatenate([np.conj(b[::-1]), np.zeros(n_pad)]))
    r_f = a_f*b_f
    r = np.fft.ifft(r_f)

    # Lags
    max_lag = len(a) - 1
    l = np.arange(-max_lag, max_lag+1)

    # Done
    return r, l

def xcorr_via_fftw(a, b):

    # Works like xcorr() in MATLAB
    # Uses FFT for convolution, via FFTW

    # Zero pad one of the inputs so they are the same length
    n_pad = len(a) - len(b)
    if n_pad > 0:
        b = np.concatenate([b, np.zeros(n_pad)])
    elif n_pad < 0:
        a = np.concatenate([a, np.zeros(-n_pad)])

    # Perform forward FFT, multiply, and inverse FFT
    # Note that we need to zero pad since multiplication of the FFTs is circular
    # convolution but we want linear convolution

    # NOTE: this is actually slower than just using numpy. If the lengths of the
    # inputs are changing we'd have a lot of bookkeeping to do to make it more
    # efficient. So using numpy in this case is probably easier.

    n_pad = len(a) - 1

    a_in = pyfftw.empty_aligned(len(a) + n_pad, dtype='complex128')
    a_out = pyfftw.empty_aligned(len(a) + n_pad, dtype='complex128')
    a_fft_obj = pyfftw.FFTW(a_in, a_out)
    a_in[:] = 0
    a_in[:len(a)] = a

    b_in = pyfftw.empty_aligned(len(b) + n_pad, dtype='complex128')
    b_out = pyfftw.empty_aligned(len(b) + n_pad, dtype='complex128')
    b_fft_obj = pyfftw.FFTW(b_in, b_out)
    b_in[:] = 0
    b_in[:len(b)] = np.conj(b[::-1])

    r_in = pyfftw.empty_aligned(len(a) + n_pad, dtype='complex128')
    r_out = pyfftw.empty_aligned(len(a) + n_pad, dtype='complex128')
    r_fft_obj = pyfftw.FFTW(r_in, r_out, direction='FFTW_BACKWARD')

    a_f = a_fft_obj()
    b_f = b_fft_obj()
    r_in[:] = a_f*b_f
    r = r_fft_obj()

    # Lags
    max_lag = len(a) - 1
    l = np.arange(-max_lag, max_lag+1)

    # Done
    return r, l

# Initial test, just to make sure things match

a = np.random.randn(900) + 1j* np.random.randn(900)
b = np.random.randn(400) + 1j* np.random.randn(400)

r1, l1 = xcorr(a, b)
r2, l2 = xcorr_via_conv(a, b)
r3, l3 = xcorr_via_fft(a, b)
r4, l4 = xcorr_via_fftw(a, b)

print(np.sum(r1-r2), np.sum(r1-r3), np.sum(r1-r4), np.sum(l1-l2), np.sum(l1-l3), np.sum(l1-l4))

# Now benchmark

# Test parameters
n_samp_a = 1250
n_samp_b = 1000
n_trials = 1000

a = np.random.randn(n_samp_a) + 1j* np.random.randn(n_samp_a)
b = np.random.randn(n_samp_b) + 1j* np.random.randn(n_samp_b)

# Basic xcorr
start_time = time.time()
for i in range(n_trials):
    r, l = xcorr(a, b)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('xcorr: %.3f ms' % ms)

# xcorr using convolution
start_time = time.time()
for i in range(n_trials):
    r, l = xcorr_via_conv(a, b)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('xcorr with conv: %.3f ms' % ms)

# xcorr using FFT
start_time = time.time()
for i in range(n_trials):
    r, l = xcorr_via_fft(a, b)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('xcorr with FFT: %.3f ms' % ms)

# xcorr using FFT via FFTW
start_time = time.time()
for i in range(n_trials):
    r, l = xcorr_via_fftw(a, b)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('xcorr with FFT from FFTW: %.3f ms' % ms)