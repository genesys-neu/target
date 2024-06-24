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

% Tests channel model

clear

% Source node parameters
src_params.id = 0;
src_params.xyz = [10, 20, 30];
src_params.fs = 1e6;
src_params.time = 3;
src_params.freq = 500e6;
src_params.freq_off = 321;
src_params.tx_phase_off = -81;
src_params.rx_phase_off = 117;
src_params.noise_fig = 10;
src_params.pow = 13;

% Destination node parameters
dst_params.id = 1;
dst_params.xyz = [110, 200, 300];
dst_params.fs = 1e6;
dst_params.time = 3.523453443113;
dst_params.freq = 500e6;
dst_params.freq_off = 321;
dst_params.tx_phase_off = -81;
dst_params.rx_phase_off = 117;
dst_params.noise_fig = 10;
dst_params.pow = 13;

% Create the nodes and channel
src = node(src_params);
dst = node(dst_params);
chan = channel(src, dst);

% Create a signal
n = 1000;
x = randn(1, n) + 1j*randn(1, n);

% Simulate the channel
y = chan.sim(x);

% Verify path length
path_len = sqrt((src_params.xyz(1) - dst_params.xyz(1))^2 + ...
                (src_params.xyz(2) - dst_params.xyz(2))^2 + ...
                (src_params.xyz(3) - dst_params.xyz(3))^2);
fprintf('Expected length : %.4f m\n', path_len);
fprintf('Actual length   : %.4f m\n', chan.path_len);

% Verify path loss
path_loss = 10*log10(var(x)/var(y));
fprintf('Expected loss   : %.4f dB\n', chan.path_loss);
fprintf('Actual loss     : %.4f dB\n', path_loss);

% Verify delay and phase
n_interp = 100;
[r, lags] = xcorr(resample(y, n_interp, 1), resample(x, n_interp, 1));
[~, idx] = max(abs(r));
phase = rad2deg(angle(r(idx)));
delay = lags(idx)/n_interp;
fprintf('Expected phase  : %.4f degrees\n', chan.phase);
fprintf('Actual phase    : %.4f degrees\n', phase);
fprintf('Expected delay  : %.4f samples\n', chan.n_delay);
fprintf('Actual delay    : %.4f samples\n', delay);