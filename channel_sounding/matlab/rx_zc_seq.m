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

function [t, mag, phase] = rx_zc_seq(filename)
%RX_ZC_SEQ Process the received Zadoff-Chu sequence

%   t        : vector of sample times (seconds)
%   mag      : estimated signal magnitude
%   phase    : estimated signal phase (radians)
%   filename : input file

% Status
fprintf('Processing %s\n', filename);

% Sample rate (SPS)
fs = 1e6;

% Refrence samples
x = read_cshort_binary('zc.dat');

% Read in the received data
fprintf('Reading data...\n');
y = read_cshort_binary(filename);
fprintf('Read in %i samples (%.4f seconds)\n', length(y), length(y)/fs);

% First find the delay in the received waveform via cross-correlation
r = xcorr(x, y(1:length(x)));
[~, max_idx] = max(abs(r).^2);
delay = length(x) - max_idx;

% Remove the delay
y_align = y(delay+1:end);

% Truncate so it's the same size as an integer number of repetitions
n_rep = floor(length(y_align)/length(x));
y_align = y_align(1:(n_rep*length(x)));

% Can split into segments and correlate each of them
% Found that this doesn't work well so leaving it set to 1 segment
% (correlate the entire code at once)
n_seg = 1;

% Cross-correlate and use the peak to find the magnitude and phase
mag = zeros(1, n_rep*n_seg);
phase = zeros(1, n_rep*n_seg);
fprintf('Running %i cross-correlations...\n', n_rep);
for i = 0:(n_rep-1)
    idx = (1:length(x)) + i*length(x);
    y_i = y_align(idx);
    for j = 0:(n_seg-1)
        idx = (1:length(x)/n_seg) + j*length(x)/n_seg;
        r = xcorr(x(idx), y_i(idx));
        [max_val, max_idx] = max(abs(r).^2);
        mag(i*n_seg+1+j) = 10*log10(max_val);
        phase(i*n_seg+1+j) = angle(r(max_idx));
        keyboard
    end
end
fprintf('Done\n');

% Unwrap the phase
phase = unwrap(phase);

% Make the time axis
n = 0:(n_rep*n_seg-1);
t_rep = length(x)/(fs/n_seg);
t = n*t_rep;

end