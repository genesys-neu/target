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

% Tests BPSK modulator

clear

% BPSK modulator parameters
alpha = 0.3;
span = 12;
sps = 8;

% Create the BPSK modulator
md = bpsk_mod(alpha, span, sps);

% Create input bits
n = 100;
bits = rand(1, n) > 0.5;

% Modulate
x = md.sim(bits);

% Verify output power
fprintf('Signal power: %.4f W\n', var(x))

% Plot
figure
plot(real(x), '-o')
grid on