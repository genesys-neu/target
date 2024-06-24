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

% Quick look at RX data

clear

%% Read in data

% Get file from browser
file_path = uigetdir('../data/');
file = fullfile(file_path, 'rx.dat');
[~, collect_time] = fileparts(file_path);
collect_time = str2double(collect_time);

% Load and parse metadata
meta_file = strcat(file, '.hdr.mat');
meta = load(meta_file);
if length(meta.rx_time) > 1
    warning('Detected dropped samples');
end
fs = meta.rx_rate(1);
fc = meta.rx_freq(1);

% Read samples
y = read_cshort_binary(file);

%% Plots

% Time axis
n = 0:(length(y)-1);
ts = 1/fs;
t = n*ts;

% Time domain
figure
plot(t, abs(y).^2)
grid on
xlabel('Time (sec)')
ylabel('Power')
title(sprintf('Time = %i, Freq = %.0f MHz, Rate = %.0f MSPS', collect_time, fc/1e6, fs/1e6))

% Frequency domain (averaged)
figure
n_fft = 1024*4;
n_overlap = 0;
win = blackmanharris(n_fft);
pwelch(y, win, n_overlap, n_fft, fs, 'centered');
title(sprintf('Time = %i, Freq = %.0f MHz, Rate = %.0f MSPS', collect_time, fc/1e6, fs/1e6))

% Spectrogram
figure
spectrogram(y, win, n_overlap, n_fft, fs, 'centered')