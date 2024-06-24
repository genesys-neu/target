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

% Checks the RX data from the debug files

clear

% Specify which node
node_id = 0;

% Load file and metadata
iq_file = sprintf('rx_msg_%i.dat', node_id);
y = read_complex_binary(iq_file);
meta_file = sprintf('rx_msg_%i.dat.hdr.mat', node_id);
meta = load(meta_file);

% Separate packets
y_pkts = {};
idx = 0;
for i = 1:length(meta.n_items)
    if meta.n_items(i) > 0
        idxs = idx + (1:meta.n_items(i));
        y_pkts{i} = y(idxs);
        idx = idxs(end);
    end
end
n_pkts = length(y_pkts);

%% Some basic plots

figure
plot(abs(y).^2)
grid on
title('Power')

figure
hold on
for i = 1:n_pkts
    plot(abs(y_pkts{i}).^2)
end
hold off
grid on
title('Power per packet')

figure
plot(y, '.')
grid on
title('Constellation')

figure
hold on
plot(real(y))
plot(imag(y))
hold off
grid on
title('Real and Imag')

%% FFT estimator for frequency offset (matches real-time code)

rrc_filt = rcosdesign(0.35, 11, 4, 'sqrt');
n_dec = 100;
fs = 1e6/n_dec;
n_fft = 2^16;
f = (fs/2)*(-1:(2/n_fft):(1-2/n_fft));
f_est = zeros(1, n_pkts);
n_avg = 1;
f_est_buff = zeros(1, n_avg);
buff_empty = 1;
buff_idx = 1;
for i = 1:n_pkts
    y_in = y_pkts{i};
    if (var(y_in) > 1e-3) && (length(y_in) >= 1000)
        
        y_in = y_in - mean(y_in);

        y_sq = y_in.^2;
        y_sq = resample(y_sq, 1, n_dec);
        y_f = fftshift(fft(y_sq, n_fft));
        [~, max_idx] = max(abs(y_f));
        
        if buff_empty
            f_est_buff(:) = f(max_idx)/2;
            buff_empty = 0;
        else
            f_est_buff(buff_idx) = f(max_idx)/2;
        end
        buff_idx = buff_idx + 1;
        if buff_idx > n_avg
            buff_idx = 1;
        end
        f_est(i) = mean(f_est_buff);
        
    else
        f_est(i) = NaN;
    end
end
fprintf('Mean %.3f, Var %.3f\n', mean(f_est), var(f_est));

figure
plot(f_est, 'o')
grid on
title('Frequency Offset')