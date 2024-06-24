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

% Processes the leader and follower logs from experiments

clear

% Path to the files
leader_file = 'leader_log.txt';
follower_2_file = 'follower_2_log.txt';
follower_3_file = 'follower_3_log.txt';

% Set the min/max run time for processing
min_t = 0;
max_t = inf;

%% Leader

% Save time and the two follower power levels
t_l = [];
p_f2 = [];
p_f3 = [];

fid = fopen(leader_file);
l = fgets(fid);
n = 1;
while 1
    
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
    
    t_l = [t_l, str2double(l{1})];
    p_f2 = [p_f2, str2double(l{6})];
    p_f3 = [p_f3, str2double(l{8})];
    
end

% Reduce to the specified times
[~, s_idx] = min(abs(t_l - min_t));
[~, e_idx] = min(abs(t_l - max_t));
idx = s_idx:e_idx;
%idx = 1:length(t_l); s_idx = 1;

%% Leader plots

figure

subplot(2,1,1)
plot(t_l(idx) - t_l(s_idx), p_f2(idx), 'o', 'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410])
grid on
legend('Follower 1', 'interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

subplot(2,1,2)
plot(t_l(idx) - t_l(s_idx), p_f3(idx), 'o', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980])
grid on
legend('Follower 2', 'location', 'se', 'interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

h = axes(gcf, 'visible', 'off');
h.XLabel.Visible='on';
h.YLabel.Visible = 'on';
xlabel(h, 'Time (seconds)', 'interpreter', 'latex')
yl = ylabel(h, 'Phase Received at Leader (degrees)', 'interpreter', 'latex');
yl.Position(1) = -0.1;
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

%% Follower 2

% Save time and frequency offset
t_f2 = [];
f_f2 = [];

fid = fopen(follower_2_file);
l = fgets(fid);
n = 1;
while 1
    
    l = fgets(fid);
    n = n + 1;
    if ~ischar(l)
        break
    end
    
    l = split(l, ', ');
    
    if length(l) < 9
        fprintf('Skipping line %i: %s\n', n, l{1});
        continue
    end
    
    t_f2 = [t_f2, str2double(l{1})];
    f_f2 = [f_f2, str2double(l{4})];
    
end

% Reduce to the specified times
[~, f2_s_idx] = min(abs(t_f2 - min_t));
[~, f2_e_idx] = min(abs(t_f2 - max_t));
f2_idx = f2_s_idx:f2_e_idx;
%f2_idx = 1:length(t_f2); f2_s_idx = 1;

%% Follower 3

% Save time and frequency offset
t_f3 = [];
f_f3 = [];

fid = fopen(follower_3_file);
l = fgets(fid);
n = 1;
while 1
    
    l = fgets(fid);
    n = n + 1;
    if ~ischar(l)
        break
    end
    
    l = split(l, ', ');
    
    if length(l) < 9
        fprintf('Skipping line %i: %s\n', n, l{1});
        continue
    end
    
    t_f3 = [t_f3, str2double(l{1})];
    f_f3 = [f_f3, str2double(l{4})];
    
end

% Reduce to the specified times
[~, f3_s_idx] = min(abs(t_f3 - min_t));
[~, f3_e_idx] = min(abs(t_f3 - max_t));
f3_idx = f3_s_idx:f3_e_idx;
%f3_idx = 1:length(t_f3); f3_s_idx = 1;

%% Follower plots

figure

subplot(2,1,1)
plot(t_f2(f2_idx) - t_f2(f2_s_idx), f_f2(f2_idx), 'o', 'Color', [0, 0.4470, 0.7410], 'MarkerFaceColor', [0, 0.4470, 0.7410])
grid on
legend('Follower 1', 'location', 'se', 'interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

subplot(2,1,2)
plot(t_f3(f3_idx) - t_f3(f3_s_idx), f_f3(f3_idx), 'o', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerFaceColor', [0.8500, 0.3250, 0.0980])
grid on
legend('Follower 2', 'interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', 16, 'FontWeight', 'bold');

h = axes(gcf, 'visible', 'off');
h.XLabel.Visible='on';
h.YLabel.Visible = 'on';
xlabel(h, 'Time (seconds)', 'interpreter', 'latex')
yl = ylabel(h, 'Frequency Offset From Leader (Hz)', 'interpreter', 'latex');
yl.Position(1) = -0.1;
set(gca, 'FontSize', 16, 'FontWeight', 'bold');