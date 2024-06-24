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

% Tests Gold code generator

clear

% Make two code generators
gc0 = gold_code_gen(0, 1000);
gc1 = gold_code_gen(1, 1000);

% Generate bits from each
b0 = gc0.gen_bits();
b1 = gc1.gen_bits();

% Map to symbols
s0 = 1 - 2*b0;
s1 = 1 - 2*b1;

% Auto and cross correlations
a0 = xcorr(s0, s0);
a1 = xcorr(s1, s1);
c01 = xcorr(s0, s1);
c10 = xcorr(s1, s0);

% Plots
figure
subplot(2,2,1)
plot(a0)
grid on
title('Autocorrelation Code 0')
subplot(2,2,2)
plot(a1)
grid on
title('Autocorrelation Code 1')
subplot(2,2,3)
plot(c01)
grid on
title('Crosscorrelation Code 0/1')
subplot(2,2,4)
plot(c10)
grid on
title('Crosscorrelation Code 1/0')