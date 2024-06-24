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

function [t, mag, phase] = rx_tone(filename)
%RX_TONE Process the received tone

%   t        : vector of sample times (seconds)
%   mag      : estimated signal magnitude
%   phase    : estimated signal phase (radians)
%   filename : input file

% Status
fprintf('Processing %s\n', filename);

% Sample rate (SPS)
fs = 1e6;

% Refrence samples
x = read_cshort_binary('tone.dat');

% Read in the received data
fprintf('Reading data...\n');
y = read_cshort_binary(filename);
fprintf('Read in %i samples (%.4f seconds)\n', length(y), length(y)/fs);

%{
% Uncomment to see the raw PSD and design the filter bandwidth. Based on
% doing this for each frequency, we pass 6-14 kHz.
fprintf('Plotting PSD...\n')
figure
pwelch(y, [], 0, 1e6, fs, 'centered')
keyboard
%}

% Bandpass filter to isolate the tone, designed using FDATOOL
Fs = 1000000;              % Sampling Frequency
Fstop1 = 5000;             % First Stopband Frequency
Fpass1 = 6000;             % First Passband Frequency
Fpass2 = 14000;            % Second Passband Frequency
Fstop2 = 15000;            % Second Stopband Frequency
Dstop1 = 0.0001;           % First Stopband Attenuation
Dpass  = 0.0057563991496;  % Passband Ripple
Dstop2 = 0.0001;           % Second Stopband Attenuation
dens   = 20;               % Density Factor

% Calculate the order from the parameters using FIRPMORD
[N, Fo, Ao, W] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2]/(Fs/2), ...
                          [0 1 0], [Dstop1 Dpass Dstop2]);

% Calculate the coefficients using the FIRPM function
b  = firpm(N, Fo, Ao, W, {dens});

% Filter the received data and remove the transients
fprintf('Filtering data...\n');
n_pad = floor((length(b)-1)/2);
y_pad = [y; zeros(n_pad, 1)];
y_filt = filter(b, 1, y_pad);
y_filt(1:n_pad) = [];

% Repeat the reference signal to match the received data
n_rep = ceil(length(y_filt)/length(x));
x_rep = repmat(x, 1, n_rep);
x_rep = x_rep(1:length(y_filt)).';

% Create a time vector
ts = 1/fs;
n = 0:(length(y_filt)-1);
t = n*ts;

% Compute magnitude of received data
mag = abs(y_filt);

% Compute phase between reference and received data
fprintf('Computing phase...\n');
phase = unwrap(angle(x_rep.*conj(y_filt)));
fprintf('Done\n');

end