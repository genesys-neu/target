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

%% Contribution of the combined phase errors

% This includes the phase estimation error terms

clear

% Number of array elements
N = [2, 3, 5, 10];

% Standard deviation of phase estimation errors (degrees)
% This applies to the leader and followers, and all signals
sigma_phase = 0:1:90;

% Threshold for normalized gain (we will look at the probability of gain
% being greater than or equal to this value)
norm_bf_gain_thresh = 0.9;

% Number of Monte Carlo trials
n_trials = 1000;

% Preallocate
n_N = length(N);
n_phases = length(sigma_phase);
prob_bf_gain = zeros(n_N, n_phases);
bf_gain = zeros(1, n_trials);

% Sweep number of array elements
for i = 1:n_N
    
    % Maximum gain for this size array
    max_bf_gain = N(i)^2;
    
    % Sweep standard deviation of phase error
    for j = 1:n_phases
        
        % Status
        fprintf('N = %i, sigma_phase = %g\n', N(i), sigma_phase(j));

        % Run Monte Carlo trials
        for k = 1:n_trials
            
            % Sample phase errors. Leader has one term, followers have
            % three terms.
            phases = [sigma_phase(j)*randn, ...
                      sum(sigma_phase(j)*randn(3, N(i) - 1))];
            
            % Beamforming gain
            bf_gain(k) = abs(sum(exp(1j*deg2rad(phases))))^2;   

        end

        % Normalize beamforming gains
        norm_bf_gain = bf_gain/max_bf_gain;
        
        % Compute the probability of beamforming gain >= our threshold
        prob_bf_gain(i, j) = ...
            sum(norm_bf_gain >= norm_bf_gain_thresh)/n_trials;

    end
    
end

% Plot one line for each array size
figure
hold on
for i = 1:n_N
    plot(sigma_phase, prob_bf_gain(i, :), 'DisplayName', num2str(N(i)))
end
hold off
grid on
xlim([min(sigma_phase), max(sigma_phase)])
ylim([0, 1])
xlabel('\sigma_{\theta} (degrees)')
ylabel(sprintf('P(G >= %g)', norm_bf_gain_thresh))
lgd = legend;
title(lgd, 'N')
boldify

%% Contribution of frequency errors

% Standard deviation of phase errors is fixed from the result above. We
% must also choose a packet time and duration.

clear

% Number of array elements
N = [2, 3, 5, 10];

% Standard deviation of frequency errors (Hz)
sigma_freq = 0:0.2:6;

% Standard deviation of phase errors (degrees)
sigma_phase = 5;

% Packet start time (seconds)
% For modeling the phase error due to frequency error
t_start = 49e-3;

% Packet duration (seconds)
t_dur = 1e-3;

% Threshold for normalized gain (we will look at the probability of gain
% being greater than or equal to this value)
norm_bf_gain_thresh = 0.9;

% Number of Monte Carlo trials
n_trials = 1000;

% Sample rate (SPS)
fs = 1e6;

% Preallocate
n_N = length(N);
n_freqs = length(sigma_freq);
prob_bf_gain = zeros(n_N, n_freqs);
bf_gain = zeros(1, n_trials);

% Create time vector
ts = 1/fs;
t = 0:ts:(t_dur - ts);

% Sweep number of array elements
for i = 1:n_N
    
    % Maximum gain for this size array
    max_bf_gain = N(i)^2;
    
    % Sweep standard deviation of frequency error
    for j = 1:n_freqs
        
        % Status
        fprintf('N = %i, sigma_freq = %g\n', N(i), sigma_freq(j));

        % Run Monte Carlo trials
        for k = 1:n_trials

            % Sample frequency errors
            % Note: leader is the true frequency (0) and the followers have
            % random offsets around it.
            freqs = [0, sigma_freq(j)*randn(1, N(i) - 1)];
            
            % Signals with frequency offsets. We can think of these as
            % being applied to a constant value (one symbol).
            sigs = exp(1j*2*pi*freqs'.*t);
                        
            % Sample phase errors. Leader has one term, followers have
            % three terms.
            phases = [sigma_phase*randn, ...
                      sum(sigma_phase*randn(3, N(i) - 1))];
                  
            % Phase errors due to frequency errors at the packet start time
            pf_terms = exp(-1j*2*pi*freqs*t_start);
            
            % Signals with phase offsets and errors
            sigs = sigs.*exp(1j*deg2rad(phases')).*(pf_terms.');
            
            % Beamforming gain. This integrates over the entire signal.
            bf_gain(k) = sum(abs(sum(sigs).^2));

        end

        % Normalize beamforming gains. We normalize by the array size and
        % the duration of the signal.
        norm_bf_gain = bf_gain/(max_bf_gain*length(t));
        
        % Compute the probability of beamforming gain >= our threshold
        prob_bf_gain(i, j) = ...
            sum(norm_bf_gain >= norm_bf_gain_thresh)/n_trials;

    end
    
end

% Plot one line for each array size
figure
hold on
for i = 1:n_N
    plot(sigma_freq, prob_bf_gain(i, :), 'DisplayName', num2str(N(i)))
end
hold off
grid on
xlim([min(sigma_freq), max(sigma_freq)])
ylim([0, 1])
xlabel('\sigma_f (Hz)')
ylabel(sprintf('P(G >= %g)', norm_bf_gain_thresh))
lgd = legend;
title(lgd, 'N')
title(sprintf('Start = %g (ms), Duration = %g (ms)', t_start*1e3, t_dur*1e3));
boldify

%% Full model, sweeping target position over grid

% Assumptions:
%   - All nodes have equal transmit powers
%   - One symbol period is >> (delay spread + timing errors), so the delay
%     on the baseband waveform is not modeled (in practice this could cause
%     ISI and require equalization)

clear

% Target grid (m)
% Altitude is zero
x_grid = -50:10:50;
y_grid = -50:10:50;

% Array element positions (m)
% Altitude is fixed for all elements
% First element is the leader
if 1
    % Line array
    p_xyz = [  0, 0, 10; ...
             -20, 0, 10; ...
             -10, 0, 10; ...
              10, 0, 10; ...
              20, 0  10];
else
    % Square array
    p_xyz = [  0,  0,  10; ...
             -10,  10, 10; ...
              10,  10, 10; ...
             -10, -10, 10; ...
              10, -10, 10];
end
      
% Center frequency (Hz)
% This is the frequency the leader uses, and followers sync to
fc = 915e6;

% Frequency mismatch between leader/target (Hz)
freq_mismatch = 100e3;

% Standard deviation of frequency errors (Hz)
sigma_freq = 0.5;

% Standard deviation of phase errors (degrees)
sigma_phase = 5;

% Packet start time (seconds)
% For modeling the phase error due to frequency error
t_start = 49e-3;

% Packet duration (seconds)
t_dur = 1e-3;

% Sample rate (SPS)
fs = 1e6;

% Number of Monte Carlo trials
n_trials = 100;
 
% Preallocate
n_x = length(x_grid);
n_y = length(y_grid);
bf_gain = zeros(1, n_trials);
avg_bf_gain = zeros(n_x, n_y);
delay_spread = zeros(n_x, n_y);

% Speed of light (m/s)
c = physconst('lightspeed');
            
% Wavelength (m)
lambda = c/fc;
      
% Number of array elements
N = size(p_xyz, 1);

% Maximum gain for this size array
max_bf_gain = N^2;

% Create time vector
ts = 1/fs;
t = 0:ts:(t_dur - ts);

% Sweep target position over grid
for i = 1:n_y
    for j = 1:n_x
        
        % Status
        fprintf('x = %i, y = %i\n', y_grid(i), x_grid(j));
        
        % Target position (m)
        t_xyz = [x_grid(j), y_grid(i), 0];
        
        % Path lengths (m)
        path_lens = pdist2(t_xyz, p_xyz);
        
        % Path delays (sec)
        path_delays = path_lens/c;
        
        % Delay spread (sec)
        delay_spread(i, j) = max(path_delays) - min(path_delays);
        
        % Path losses, normalized to the first element
        path_losses = 1./((4*pi*path_lens/lambda).^2);
        norm_path_losses = path_losses/path_losses(1);
        
        % Run Monte Carlo trials
        for k = 1:n_trials
            
            % Sample frequency errors
            % Note: leader is the true frequency (0) and the followers have
            % random offsets around it.
            freqs = [0, sigma_freq*randn(1, N - 1)];
            
            % Signals with frequency offsets. We can think of these as
            % being applied to a constant value (one symbol).
            sigs = exp(1j*2*pi*freqs'.*t);
            
            % Sample phase errors. Leader has one term, followers have
            % three terms.
            phases = [sigma_phase*randn, ...
                      sum(sigma_phase*randn(3, N - 1))];
            
            % Frequency mismatch between leader/target
            fm_terms = exp(1j*2*pi*(freq_mismatch + freqs).*path_delays);
            
            % Phase error terms (all combined)
            pe_terms = exp(1j*deg2rad(phases));
            
            % Phase errors due to frequency errors at the packet start time
            pf_terms = exp(-1j*2*pi*freqs*t_start);
            
            % Frequency errors between leader/follower
            f_terms = exp(1j*2*pi*freqs.*path_delays);
            
            % Combine terms
            comb_terms = norm_path_losses'.*...
                         sigs.*...
                         (fm_terms.').*...
                         (f_terms.').*...
                         (pe_terms.').*...
                         (pf_terms.');
            
            % Beamforming gain. This integrates over the entire signal.
            bf_gain(k) = sum(abs(sum(comb_terms)).^2);
            
        end
        
        % Single element gain, integrating over the entire signal
        se_gain = sum(abs(comb_terms(1, :)).^2);
        
        % Average gain by doing beamforming
        avg_bf_gain(i, j) = mean(bf_gain)/se_gain;
        
    end
end
%% Plots

% Mean gain over the grid
figure
imagesc(x_grid, y_grid, 10*log10(avg_bf_gain))
hold on
for i = 1:N
    plot(p_xyz(i, 1), p_xyz(i, 2), 'sw', 'MarkerFaceColor', 'w')
end
hold off
set(gca, 'YDir', 'normal')
cbar = colorbar;
xlabel('x (meters)')
ylabel('y (meters)')
ylabel(cbar, 'Average Gain (dB)', 'FontSize', 16, 'FontWeight', 'bold')
boldify
% Inner dot without boldify
hold on
for i = 1:N
    plot(p_xyz(i, 1), p_xyz(i, 2), 'sw', 'MarkerFaceColor', 'k')
end
hold off

% Delay spread over grid
figure
imagesc(x_grid, y_grid, delay_spread*1e6)
hold on
for i = 1:N
    plot(p_xyz(i, 1), p_xyz(i, 2), 'sw', 'MarkerFaceColor', 'w')
end
hold off
set(gca, 'YDir', 'normal')
cbar = colorbar;
xlabel('x (meters)')
ylabel('y (meters)')
ylabel(cbar, 'Delay Spread (usec)', 'FontSize', 16, 'FontWeight', 'bold')
boldify
% Inner dot without boldify
hold on
for i = 1:N
    plot(p_xyz(i, 1), p_xyz(i, 2), 'sw', 'MarkerFaceColor', 'k')
end
hold off