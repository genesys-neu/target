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

% Generate a tone (complex sinusoid) for narrowband testing

clear

% Sample rate (SPS)
fs = 1e6;

% Baseband tone frequency (Hz)
fc = 10e3;

% One period of samples
t_dur = 1/fc;
ts = 1/fs;
t = 0:ts:(t_dur - ts);
x = 0.9*exp(1j*2*pi*fc*t);

% Save to file
write_cshort_binary(x, 'tone.dat');

% Plot
figure
hold on
plot(real(x), '-o')
plot(imag(x), '-o')
hold off
grid on
xlabel('Samples')
ylabel('Amplitude')
legend('Real', 'Imag')
title('One period of complex sinusoid')

% Next, make sure there's no phase glitch when the SDR plays it out

% Repeat the signal
n_rep = 1000;
x_rep = repmat(x, 1, n_rep);

% Make an independent copy based on the longer time duration, and give it a
% lower amplitude and phase shift
t_long_dur = n_rep*(1/fc);
t_long = 0:ts:(t_long_dur - ts);
phase_offset = pi/7;
y = 0.1*exp(1j*2*pi*fc*t_long)*exp(1j*phase_offset);

% Compute the phase difference
phase_diff = diff(angle(x_rep.*conj(y)));

% Plot
figure
plot(rad2deg(phase_diff), '-o')
grid on
xlabel('Samples')
ylabel('Phase Difference (degrees)')