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

# Benchmarks FFT to see which is the fastest (numpy, scipy, FFTW)

import pyfftw
import numpy as np
import time
import scipy as sp

# Test parameters
n_samp = 1000
n_fft = 2**16
n_trials = 1000

# Make data
x = np.random.randn(n_samp) + 1j*np.random.randn(n_samp)

# Basic numpy
start_time = time.time()
for i in range(n_trials):
	y = np.fft.fft(x, n_fft)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('numpy: %.3f ms' % ms)

# numpy using FFTW interface
x_in = pyfftw.empty_aligned(n_samp, dtype='complex128')
pyfftw.interfaces.cache.enable()
x_in[:] = x
start_time = time.time()
for i in range(n_trials):
	y = pyfftw.interfaces.numpy_fft.fft(x_in, n_fft)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('FFTW via numpy: %.3f ms' % ms)

# Basic scipy
start_time = time.time()
for i in range(n_trials):
	y = sp.fft.fft(x, n_fft)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('scipy: %.3f ms' % ms)

# scipy using FFTW interface
x_in = pyfftw.empty_aligned(n_samp, dtype='complex128')
pyfftw.interfaces.cache.enable()
x_in[:] = x
start_time = time.time()
for i in range(n_trials):
	y = pyfftw.interfaces.scipy_fft.fft(x_in, n_fft)
end_time = time.time()
ms = (end_time - start_time)*1e3
print('FFTW via scipy: %.3f ms' % ms)

# FFTW on its own
x_in = pyfftw.empty_aligned(n_fft, dtype='complex128')
y_out = pyfftw.empty_aligned(n_fft, dtype='complex128')
x_in[:n_samp] = x
fft_obj = pyfftw.FFTW(x_in, y_out)
start_time = time.time()
for i in range(n_trials):
	y = fft_obj()
end_time = time.time()
ms = (end_time - start_time)*1e3
print('FFTW: %.3f ms' % ms)