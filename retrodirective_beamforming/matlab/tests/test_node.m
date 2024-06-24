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

% Tests node object

clear

% Node parameters
params.id = 0;
params.xyz = [10, 20, 30];
params.fs = 1e6;
params.time = 123.43;
params.freq = 500e6;
params.freq_off = 321;
params.tx_phase_off = -81;
params.rx_phase_off = 117;
params.noise_fig = 10;
params.pow = 13;

% Create the node
m = node(params);

% Make sure we can save data
m.data.a = 1;
m.data.b = 2;
m.data.c = 3;

% Display the node properties
disp(m);
disp(m.tx);
disp(m.rx);
disp(m.data);