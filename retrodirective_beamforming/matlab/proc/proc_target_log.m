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

% Processes the taget log from experiments

clear

% Path to the file
file_in = 'target_log.txt';

% Set the min/max run time for processing
min_t = 0;
max_t = inf;

%% Target

t = [];
t_l = [];
l_pow = [];
t_f = [];
f_pow = [];
t_bf = [];
bf_pow = [];

fid = fopen(file_in);
l = fgets(fid);
l = split(l, ', ');
l_prev1 = l;
n = 1;
while 1
    
    l_prev2 = l_prev1;
    l_prev1 = l;
    l = fgets(fid);
    n = n + 1;
    if ~ischar(l)
        break
    end
    
    l = split(l, ', ');
    
    if length(l) == 1
        fprintf('Skipping line %i: %s\n', n, l{1});
        continue
    end
    
    if strcmp(l{1}(1:3), 'exp')    
        if length(l_prev1) > 1 && length(l_prev2) > 1
            if strcmp(l_prev1{3}(1:3), 'lea') && strcmp(l_prev2{3}(1:3), 'fol')
                
                t_l = [t_l, str2double(l_prev1{1})];
                l_pow = [l_pow, str2double(l_prev1{2})];
                
                t_f = [t_f, str2double(l_prev2{1})];
                f_pow = [f_pow, str2double(l_prev2{2})];
                
                t_bf = [t_bf, str2double(l_prev1{1}) + 0.050];
                bf_pow = [bf_pow, str2double(l{2}(5:end))];
                
            end
        end
    end
    
end

% Calculate the optimal beamformed power level
opt_bf_pow = (sqrt(l_pow) + sqrt(f_pow)).^2;

[~, s_idx] = min(abs(t_bf - min_t));
[~, e_idx] = min(abs(t_bf - max_t));
idx = s_idx:e_idx;
%idx = 1:length(t_bf); s_idx = 1;

%% Results

% Plot of all power levels
figure
hold on
plot(t_f(idx) - t_f(s_idx), 10*log10(f_pow(idx)), 'o', 'MarkerFaceColor', [0, 0.4470, 0.7410])
plot(t_l(idx) - t_l(s_idx), 10*log10(l_pow(idx)), 'o', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980])
plot(t_bf(idx) - t_bf(s_idx), 10*log10(bf_pow(idx)), 'o', 'MarkerFaceColor', [0.9290, 0.6940, 0.1250])
plot(t_bf(idx) - t_bf(s_idx), 10*log10(opt_bf_pow(idx)), 'k', 'LineWidth', 3)
hold off
grid on
xlabel('Time (seconds)', 'interpreter', 'latex');
ylabel('Received Power (dBFS)', 'interpreter', 'latex');
legend('Follower', 'Leader', 'Beamforming', 'Optimal Beamforming', 'location', 'se', 'interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

% Percent of optimal power level achieved
p_bf = bf_pow./opt_bf_pow;

% Probability of beamformed power level >= 90%
sum(p_bf(idx)>=0.9)/length(p_bf(idx))

% Mean of beamformed power levels
mean(p_bf(idx))

% Plot percent of optimal power level over time
figure
plot(t_bf(idx) - t_bf(s_idx), p_bf(idx), 'o')
grid on