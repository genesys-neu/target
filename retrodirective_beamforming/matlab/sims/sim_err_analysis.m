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

%% Phase errors

% This assumes target phase errors are half of the inter-array phase errors

clear

% Number of array elements
N = [2, 4, 8, 16];

% Standard deviation of phase estimation errors (degrees)
sigma_phase = 0:0.2:30;

% Threshold for normalized gain (we will look at the probability of gain
% being greater than or equal to this value)
norm_bf_gain_thresh = 0.9;

% Number of Monte Carlo trials
n_trials = 100000;

% Preallocate
n_N = length(N);
n_phases = length(sigma_phase);
bf_gain = zeros(1, n_trials);
prob_bf_gain = zeros(n_N, n_phases);
avg_bf_gain = zeros(n_N, n_phases);

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
            
            % Leader term
            l_phase = sigma_phase(j)*randn;
            
            % Follower terms
            f_phase = sigma_phase(j)*randn(1, N(i) - 1) + ...
                      sqrt(2)*sigma_phase(j)*randn(1, N(i) - 1) + ...
                      sqrt(2)*sigma_phase(j)*randn(1, N(i) - 1);
                  
            % Beamforming gain
            bf_gain(k) = abs(exp(1j*deg2rad(l_phase)) + ...
                             sum(exp(1j*deg2rad(f_phase))))^2;
            
        end

        % Normalize beamforming gains
        norm_bf_gain = bf_gain/max_bf_gain;
        
        % Compute the probability of beamforming gain >= our threshold
        prob_bf_gain(i, j) = ...
            sum(norm_bf_gain >= norm_bf_gain_thresh)/n_trials;
        
        % Compute average beamforming gain
        avg_bf_gain(i, j) = mean(norm_bf_gain);

    end
    
end

% Save results
save('phase_bf_gain.mat')

% Plot probability of beamforming gain >= our threshold
figure
hold on
for i = 1:n_N
    plot(sigma_phase, prob_bf_gain(i, :), 'DisplayName', num2str(N(i)))
end
hold off
grid on
xlim([min(sigma_phase), max(sigma_phase)])
ylim([0, 1])
xlabel('$\sigma_{\theta}$ (degrees)', 'interpreter', 'latex')
ylabel(sprintf('$P(G \\ge %g)$', norm_bf_gain_thresh), 'interpreter', 'latex')
title(sprintf('%i trials', n_trials), 'interpreter', 'latex')
lgd = legend('location', 'ne', 'interpreter', 'latex');
title(lgd, 'N', 'interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex')
boldify
savefig('phase_prob_bf_gain.fig')

% Plot average beamforming gain
figure
hold on
for i = 1:n_N
    plot(sigma_phase, avg_bf_gain(i, :), 'DisplayName', num2str(N(i)))
end
hold off
grid on
xlim([min(sigma_phase), max(sigma_phase)])
ylim([0, 1])
xlabel('$\sigma_{\theta}$ (degrees)', 'interpreter', 'latex')
ylabel('$E[G]$', 'interpreter', 'latex')
title(sprintf('%i trials', n_trials), 'interpreter', 'latex')
lgd = legend('location', 'ne', 'interpreter', 'latex');
title(lgd, 'N', 'interpreter', 'latex')
set(gca, 'TickLabelInterpreter', 'latex')
boldify
savefig('phase_avg_bf_gain.fig')

%% Frequency errors

% Ignores the frequency modulation on the packet, only considers the phase
% offsets due to protocol delay and the random phase errors based on the
% above simulation.

clear

% Number of array elements
N = [2, 4, 8, 16];

% Standard deviation of phase estimation errors (degrees)
sigma_phase = 5;

% Standard deviation of frequency estimation errors (Hz)
sigma_freq = [0:0.1:1, 1.25:0.25:10];

% Protocol delays (sec)
delta_t = logspace(-3, 0, 40);

% Threshold for normalized gain (we will look at the probability of gain
% being greater than or equal to this value)
norm_bf_gain_thresh = 0.9;

% Number of Monte Carlo trials
n_trials = 100000;

% Preallocate
n_N = length(N);
n_freqs = length(sigma_freq);
n_times = length(delta_t);
bf_gain = zeros(1, n_trials);
prob_bf_gain = zeros(n_N, n_freqs, n_times);
avg_bf_gain = zeros(n_N, n_freqs, n_times);

% Sweep number of array elements
for i = 1:n_N
    
    % Maximum gain for this size array
    max_bf_gain = N(i)^2;
    
    % Sweep standard deviation of frequency error
    for j = 1:n_freqs
        
        % Sweep protocol delay
        for k = 1:n_times
        
            % Status
            fprintf('N = %i, sigma_freq = %g, delta_t = %g\n', N(i), sigma_freq(j), delta_t(k));

            % Run Monte Carlo trials
            for l = 1:n_trials

                % Leader term
                l_phase = sigma_phase*randn;

                % Follower terms
                f_phase = rad2deg(2*pi*sigma_freq(j)*randn(1, N(i) - 1)*delta_t(k)) + ...
                          sigma_phase*randn(1, N(i) - 1) + ...
                          sqrt(2)*sigma_phase*randn(1, N(i) - 1) + ...
                          sqrt(2)*sigma_phase*randn(1, N(i) - 1);

                % Beamforming gain
                bf_gain(l) = abs(exp(1j*deg2rad(l_phase)) + ...
                                 sum(exp(1j*deg2rad(f_phase))))^2;

            end

            % Normalize beamforming gains
            norm_bf_gain = bf_gain/max_bf_gain;

            % Compute the probability of beamforming gain >= our threshold
            prob_bf_gain(i, j, k) = ...
                sum(norm_bf_gain >= norm_bf_gain_thresh)/n_trials;

            % Compute average beamforming gain
            avg_bf_gain(i, j, k) = mean(norm_bf_gain);
            
        end

    end
    
end

% Save results
save('freq_time_bf_gain.mat')

% Plot probability of beamforming gain >= our threshold
for i = 1:n_N
    figure
    contourf(sigma_freq, delta_t, squeeze(prob_bf_gain(i, :, :))', 100, 'edgecolor', 'none')
    set(gca, 'yscale', 'log')
    xlabel('$\sigma_f$ (Hz)', 'interpreter', 'latex')
    ylabel('$T_E$ (seconds)', 'interpreter', 'latex')
    title(sprintf('N = %i, %i trials', N(i), n_trials), 'interpreter', 'latex')
    colormap(parula(10))
    cbar = colorbar;
    ylabel(cbar, sprintf('$P(G \\ge %g)$', norm_bf_gain_thresh), 'interpreter', 'latex', 'FontSize', 16)
    set(gca, 'TickLabelInterpreter', 'latex')
    set(cbar, 'TickLabelInterpreter', 'latex')
    caxis([0, 1])
    boldify
    savefig(sprintf('freq_time_prob_bf_gain_%i.fig', N(i)))
end

% Plot average beamforming gain
for i = 1:n_N
    figure
    contourf(sigma_freq, delta_t, squeeze(avg_bf_gain(i, :, :))', 'edgecolor', 'none')
    set(gca, 'yscale', 'log')
    xlabel('$\sigma_f$ (Hz)', 'interpreter', 'latex')
    ylabel('$T_E$ (seconds)', 'interpreter', 'latex')
    title(sprintf('N = %i, %i trials', N(i), n_trials), 'interpreter', 'latex')
    colormap(parula(10))
    cbar = colorbar;
    ylabel(cbar, '$E[G]$', 'interpreter', 'latex', 'FontSize', 16)
    set(gca, 'TickLabelInterpreter', 'latex')
    set(cbar, 'TickLabelInterpreter', 'latex')
    caxis([0, 1])
    boldify
    savefig(sprintf('freq_time_avg_bf_gain_%i.fig', N(i)))
end