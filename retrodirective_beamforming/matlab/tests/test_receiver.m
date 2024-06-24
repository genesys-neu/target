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

% Tests receiver object

clear

% Receiver parameters
fs = 1e6;
noise_fig = 10;
freq = 1.2e9;
freq_off = 166;
phase_off = -89;

% Create the receiver
rx = receiver(fs, noise_fig, freq, freq_off, phase_off);

% Create an input signal
ts = 1/fs;
t = 0:ts:(10e-3 - ts);
fc = 10e3;
x = 5e-6*exp(1j*2*pi*fc*t);

% Simulate the receiver at different times
t = [0, 1.33e-7, 6.3569991];
for i = 1:length(t)

    if i > 1
        fprintf('\n');
    end
    fprintf('Simulating at time %g seconds\n', t(i));
    y = rx.sim(x, t(i));

    % Signal power and SNR
    fprintf('Signal power       : %.4f dBm\n', rx.sig_pow)
    fprintf('SNR                : %.4f dB\n', rx.snr)

    % Verify noise figure
    ktb = physconst('boltzmann')*290*fs;
    act_noise_fig = 10*log10(var(rx.noise)) - 10*log10(ktb);
    fprintf('Expected noise fig : %.4f dB\n', noise_fig)
    fprintf('Actual noise fig   : %.4f dB\n', act_noise_fig);

    % Verify phase (from first sample only)
    new_phase = phase_off - rad2deg(t(i)*(freq + freq_off)*2*pi);
    new_phase = mod(new_phase, 360);
    act_phase = rad2deg(angle(y(1).*conj(x(1))));
    fprintf('Expected phase     : %.4f degrees\n', new_phase)
    fprintf('Actual phase       : %.4f degrees\n', act_phase);

    % Verify frequency error (from derivative of unwrapped phase)
    act_freq = mean(diff(unwrap(angle(y.*conj(x(:)))))*fs/(2*pi));
    fprintf('Expected freq      : %.4f Hz\n', -freq_off)
    fprintf('Actual freq        : %.4f Hz\n', act_freq);
    
end