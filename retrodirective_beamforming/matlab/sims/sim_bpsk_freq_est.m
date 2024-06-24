%*********************************************************************************
% DISTRIBUTION STATEMENT A. Approved for public release. Distribution is unlimited.
%
% This material is based upon work supported under Air Force Contract No. FA8702-15-D-0001. Any opinions, findings, conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the U.S. Air Force.
%
% (c) 2023 Massachusetts Institute of Technology.
%
% Subject to FAR52.227-11 Patent Rights - Ownership by the contractor (May 2014)
%
% The software/firmware is provided to you on an As-Is basis
%
% Delivered to the U.S. Government with Unlimited Rights, as defined in DFARS Part 252.227-7013 or 7014 (Feb 2014). Notwithstanding any copyright notice, U.S. Government rights in this work are defined by DFARS 252.227-7013 or DFARS 252.227-7014 as detailed above. Use of this work other than as specifically authorized by the U.S. Government may violate any copyrights that exist in this work.
%********************************************************************************/

% Simulates frequency offset estimation for a BPSK signal

clear

% Parameters
n_samps = 1000;
n_pad = 100;
samps_per_sym = 4;
fs = 1e6;

% Make RRC filtered BPSK
n_syms = floor(n_samps/samps_per_sym);
bits = rand(1, n_syms) > 0.5;
syms = 1 - 2*bits;
syms_up = upsample(syms, samps_per_sym);
syms_up = [syms_up, zeros(1, n_pad)];
rrc_filt = rcosdesign(0.35, 11, samps_per_sym);
x = filter(rrc_filt, 1, syms_up);

% Apply a frequency offset
freq_off = 123.456;
n = 0:(length(x) - 1);
s = exp(1j*2*pi*(freq_off/fs)*n);
x_freq_off = s.*x;

% Apply a phase offset
phase_off = deg2rad(174.36);
x_phase_off = x_freq_off*exp(1j*phase_off);

% Apply a delay
n_delay = 123;
x_delay = [zeros(1, n_delay), x_phase_off];

% Add noise
snr_dB = 300;
snr = 10^(snr_dB/10);
sig_pow = var(x);
noise_pow = sig_pow/snr;
n = length(x_delay);
noise = sqrt(noise_pow)*(randn(1, n) + 1j*randn(1, n));
y = x_delay + noise;

% Estimate frequency offset
n_fft = 2^23;
y_sq = y.^2;
y_fft = fftshift(fft(y_sq, n_fft));
f = (fs/2)*(-1:(2/n_fft):(1-2/n_fft));
[~, max_idx] = max(abs(y_fft));
freq_est = f(max_idx)/2;
fprintf('f_est = %.4f Hz\n', freq_est);

% Correct for frequency offset
y_freq = y.*exp(-1j*2*pi*(freq_est/fs)*(0:(n-1)));

% Estimate time delay and phase
[r, l] = xcorr(y_freq, x);
[~, max_idx] = max(abs(r).^2);
time_est = l(max_idx)*(1/fs);
phase_est = rad2deg(angle(r(max_idx)) + 2*pi*freq_off*time_est);
fprintf('t_est = %.4f usec\n', time_est*1e6);
fprintf('p_est = %.4f deg\n', phase_est);