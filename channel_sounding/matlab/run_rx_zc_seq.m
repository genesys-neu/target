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

% Demonstrates running the rx_zc_seq function

clear

% Get file from browser
file_path = uigetdir('../data/');
file = fullfile(file_path, 'rx.dat');

% Run the function on the file
[t, mag, phase] = rx_zc_seq(file);

% Plot the results
figure

% Magnitude
subplot(2,1,1)
plot(t, 20*log10(mag))
xlim([t(1), t(end)])
grid on
xlabel('Time (sec)')
ylabel('Magnitude (dB)')

% Phase
subplot(2,1,2)
plot(t, rad2deg(phase));
xlim([t(1), t(end)])
grid on
xlabel('Time (sec)')
ylabel('Phase (deg)')