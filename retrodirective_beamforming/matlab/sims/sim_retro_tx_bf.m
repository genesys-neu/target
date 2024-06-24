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

%% Simulation of distributed retrodirective beamforming
% Includes the basic protocol and signal processing. This is based on a
% fixed placement of the nodes, and may need modifications for arbitrary
% placements (to handle the variable length signals due to the different
% channel delays).

%% Setup

clear
clc

% Set up environment if needed
try
    tmp = gold_code_gen(0, 1);
    clear tmp
catch ME
    fprintf('Setting up environment\n\n');
    setup_env
end

%% Define parameters

% Common
params.fs = 1e6;                % Sample rate (SPS)
params.freq = 916.698451e6;     % Carrier frequency (Hz)
params.noise_fig = 10;          % Receiver noise figure (dB)
params.pow = 10;                % Transmitter power (dBm)
params.gc_len = 1000;           % Sounding gold code length (samples)
params.alpha = 0.3;             % RRC filter roll-off
params.span = 10;               % RRC filter span (symbols)
params.sps = 8;                 % RRC filter samples per symbol
params.cfc_freq_res = 1;        % Frequency estimator resolution (Hz)
params.tpe_interp = 100;        % Time/phase estimator interpolation
params.proc_delay_len = 1000;   % Processing delay (samples)
params.target_t_delay = 9.6e-6;	% Delay before target transmits (seconds)
params.rng_seed = 100;          % Random number generator seed

% Leader
params.leader.id = 0;                   % Identifier
params.leader.xyz = [0, 100, 100];      % Cartesian coordinates (meters)
params.leader.fs = params.fs;
params.leader.time = 54.648915529099;   % Clock time (seconds)
params.leader.freq = params.freq;
params.leader.freq_off = 112.698;       % Frequency offset (Hz)
params.leader.tx_phase_off = 245.258;   % Transmitter phase offset (degrees)
params.leader.rx_phase_off = -90.157;   % Receiver phase offset (degrees)
params.leader.noise_fig = params.noise_fig;
params.leader.pow = params.pow;

% Follower 1
params.followers(1).id = 1;
params.followers(1).xyz = [100, -20, 80];
params.followers(1).fs = params.fs;
params.followers(1).time = 988.248004586874;
params.followers(1).freq = params.freq;
params.followers(1).freq_off = -70.111;
params.followers(1).tx_phase_off = -181.301;
params.followers(1).rx_phase_off = 20.094;
params.followers(1).noise_fig = params.noise_fig;
params.followers(1).pow = params.pow;

% Follower 2
params.followers(2).id = 2;
params.followers(2).xyz = [200, 70, 50];
params.followers(2).fs = params.fs;
params.followers(2).time = 111.111555668010;
params.followers(2).freq = params.freq;
params.followers(2).freq_off = 55.974;
params.followers(2).tx_phase_off = 108.374;
params.followers(2).rx_phase_off = -4.344;
params.followers(2).noise_fig = params.noise_fig;
params.followers(2).pow = params.pow;

% Target
params.target.id = 3;
params.target.xyz = [0, 0, 0];
params.target.fs = params.fs;
params.target.time = 8.777000001874;
params.target.freq = params.freq;
params.target.freq_off = 201.485;
params.target.tx_phase_off = 77.006;
params.target.rx_phase_off = 154.981;
params.target.noise_fig = params.noise_fig;
params.target.pow = params.pow;

%% Initialize

% Status
fprintf('Initializing\n');

% Set random seed
rng(params.rng_seed);

% Create nodes
leader = node(params.leader);
n_followers = length(params.followers);
for i = 1:n_followers
    followers(i) = node(params.followers(i));
end
target = node(params.target);

% Create channels between nodes
for i = 1:n_followers
    l2f_chans(i) = channel(leader, followers(i));
end
t2l_chan = channel(target, leader);
for i = 1:n_followers
    t2f_chans(i) = channel(target, followers(i));
end

%{
% Plot node geometry
figure
hold on
plot3(leader.xyz(1), leader.xyz(2), leader.xyz(3), 'or')
for i = 1:n_followers
    plot3(followers(i).xyz(1), ...
          followers(i).xyz(2), ...
          followers(i).xyz(3), 'ob')
end
plot3(target.xyz(1), target.xyz(2), target.xyz(3), 'xk')
hold off
grid on
%}

%% Leader transmits to followers

% Status
fprintf('\nLeader transmits to followers\n');

% Generate and modulate bits
gc_gen = gold_code_gen(params.leader.id, params.gc_len);
bb_mod = bpsk_mod(params.alpha, params.span, params.sps);
bits = gc_gen.gen_bits();
ref_samps = bb_mod.sim(bits);

% TX simulation
% The signal starts at the leader's current time
tx_samps = leader.tx.sim(ref_samps, leader.time);

% Save data
leader.data.l2f.bits = bits;
leader.data.l2f.ref_samps = ref_samps;
leader.data.l2f.tx_samps = tx_samps;

%% Followers receive from leader

for i = 1:n_followers
    
    % Status
    if i == 1
        fprintf('\n');
    end
    fprintf('Follower %i receives from leader\n', i);
    
    % Input data
    tx_samps = leader.data.l2f.tx_samps;
    
    % Channel simulation
    chan_samps = l2f_chans(i).sim(tx_samps);
    
    % RX simulation
    % NOTE: signal starts at the follower's current time, since the output
    % of the channel includes the fractional sample delay (so the receiver
    % phase will be set correctly by the time the signal arrives).
    rx_samps = followers(i).rx.sim(chan_samps, followers(i).time);
    
    % Save data
    followers(i).data.l2f.chan_samps = chan_samps;
    followers(i).data.l2f.rx_samps = rx_samps;
    
end

%% Folowers estimate channels from leader

for i = 1:n_followers
    
    % Status
    fprintf('\nFollower %i estimates channel from leader\n', i);
    
    % Input data
    rx_samps = followers(i).data.l2f.rx_samps;
    
    % Carrier frequency offset estimation and correction
    cfc = comm.CoarseFrequencyCompensator(...
        'Modulation', 'BPSK', ...
        'FrequencyResolution', params.cfc_freq_res, ...
        'SampleRate', params.fs);
    [cfc_samps, cfc_est] = cfc(rx_samps);

    % Correlate to estimate time and phase offsets
    ref_samps = leader.data.l2f.ref_samps;
    cfc_samps_interp = resample(cfc_samps, params.tpe_interp, 1);
    ref_samps_interp = resample(ref_samps, params.tpe_interp, 1);
    [r, lags] = xcorr(cfc_samps_interp, ref_samps_interp);
    [v, idx] = max(abs(r).^2);
    n_off_est = lags(idx)/params.tpe_interp;
    phase_est = rad2deg(angle(r(idx)));
    
    % Align and estimate SNR
    tpe_samps = cfc_samps_interp(lags(idx)+1:end);
    tpe_samps = resample(tpe_samps, 1, params.tpe_interp);
    idxs = 1:length(ref_samps);
    snr_samps = tpe_samps(idxs)*exp(-1j*deg2rad(phase_est));
    snr_samps = snr_samps/std(snr_samps);
    noise_samps = snr_samps - ref_samps;
    snr_est = 10*log10(1/var(noise_samps));
    
    % Compute errors vs. expected
    cfc_est_err = cfc_est - (leader.freq_off - followers(i).freq_off);
    n_off_est_err = n_off_est - l2f_chans(i).n_delay;
    exp_phase = leader.tx_phase_off + ...
                rad2deg(leader.time*leader.tx.act_freq*2*pi) + ...
                l2f_chans(i).phase + ...
                followers(i).rx_phase_off - ...
                rad2deg(followers(i).time*followers(i).rx.act_freq*2*pi);
    exp_phase = mod(exp_phase, 360);
    phase_est_err = phase_est - exp_phase;
    snr_est_err = snr_est - ...
        (leader.pow - l2f_chans(i).path_loss - followers(i).rx.noise_pow);
    
    % Status
    fprintf('Frequency offset : %.4f Hz\n', cfc_est)
    fprintf('Frequency error  : %.4f Hz\n', cfc_est_err);
    fprintf('Time offset      : %.4f samples\n', n_off_est);
    fprintf('Time error       : %.4f samples\n', n_off_est_err);
    fprintf('Phase offset     : %.4f degrees\n', phase_est);
    fprintf('Phase error      : %.4f degrees\n', phase_est_err);
    fprintf('SNR              : %.4f dB\n', snr_est);
    fprintf('SNR error        : %.4f dB\n', snr_est_err);
    
    % Save data
    followers(i).data.l2f.cfc_est = cfc_est;
    followers(i).data.l2f.cfc_est_err = cfc_est_err;
    followers(i).data.l2f.n_off_est = n_off_est;
    followers(i).data.l2f.n_off_est_err = n_off_est_err;
    followers(i).data.l2f.phase_est = phase_est;
    followers(i).data.l2f.phase_est_err = phase_est_err;
    followers(i).data.l2f.snr_est = snr_est;
    followers(i).data.l2f.snr_est_err = snr_est_err;

end

%% All nodes delay and update clocks

% Delay includes the waveform duration as well as any protocol delay or
% processing times
n_delay = params.gc_len + params.proc_delay_len;
t_delay = n_delay/params.fs;

% Status
fprintf('\nProcessing delay, all nodes update clocks by %g seconds\n', ...
        t_delay);

% Update time for all nodes. The leader is the base time and each follower
% has an offset from its time of arrival estimate. The time is set to when
% the followers will next transmit.
leader.time = leader.time + t_delay;
for i = 1:n_followers
    n_off_est = followers(i).data.l2f.n_off_est;
    followers(i).time = followers(i).time + n_off_est/params.fs + t_delay;
end

%% Followers transmit to leader

for i = 1:n_followers
    
    % Status
    if i == 1
        fprintf('\n');
    end
    fprintf('Follower %i transmits to leader\n', i);
    
    % Generate and modulate bits
    gc_gen = gold_code_gen(params.followers(i).id, params.gc_len);
    bb_mod = bpsk_mod(params.alpha, params.span, params.sps);
    bits = gc_gen.gen_bits();
    ref_samps = bb_mod.sim(bits);
    
    % Apply estimated frequency offset
    cfc_est = followers(i).data.l2f.cfc_est;
    n = 0:(length(ref_samps) - 1);
    ref_samps_freq = ref_samps.*exp(1j*2*pi*(cfc_est/params.fs)*n');
    
    % Apply phase offset based on the protocol delay (since the time we
    % estimated the frequency offset from the leader)
    ref_samps_freq = ref_samps_freq*exp(1j*2*pi*cfc_est*2e-3);

    % TX simulation
    % The signal starts at the followers's current time
    tx_samps = followers(i).tx.sim(ref_samps_freq, followers(i).time);

    % Save data
    followers(i).data.f2l.bits = bits;
    followers(i).data.f2l.ref_samps = ref_samps;
    followers(i).data.f2l.tx_samps = tx_samps;
    
end

%% Leader receives from followers

for i = 1:n_followers
    
    % Status
    if i == 1
        fprintf('\n');
    end
    fprintf('Leader receives from follower %i\n', i);
    
    % Input data
    tx_samps = followers(i).data.f2l.tx_samps;
    
    % Channel simulation (save each signal independently)
    if i == 1
        chan_samps = l2f_chans(i).sim(tx_samps);
    else
        chan_samps = [chan_samps, l2f_chans(i).sim(tx_samps)];
    end
    
end

% Leader will receive summation of all follower signals
sum_chan_samps = sum(chan_samps, 2);

% RX simulation
% NOTE: signal starts at the leader's current time, since the output of the
% channel includes the fractional sample delay (so the receiver phase will
% be set correctly by the time the signal arrives).
rx_samps = leader.rx.sim(sum_chan_samps, leader.time);

% Save data
leader.data.f2l.chan_samps = chan_samps;
leader.data.f2l.rx_samps = rx_samps;

%% Leader estimates channels from each follower

% Input data
rx_samps = leader.data.f2l.rx_samps;

% NOTE: there is no carrier frequency offset estimation since we are
% receiving the superposition of modulated signals from the followers (we
% would need to do individual frequency estimates using the matched
% filter). Assuming the followers are now close enough in frequency so this
% doesn't matter (may need to be verified depending on the quality of the
% hardware being used).

% Interpolate for timing/phase estimation
rx_samps_interp = resample(rx_samps, params.tpe_interp, 1);

for i = 1:n_followers
    
    % Status
    fprintf('\nLeader estimates channel from follower %i\n', i);

    % Correlate to estimate time and phase offsets
    ref_samps = followers(i).data.f2l.ref_samps;
    ref_samps_interp = resample(ref_samps, params.tpe_interp, 1);
    [r, lags] = xcorr(rx_samps_interp, ref_samps_interp);
    [v, idx] = max(abs(r).^2);
    n_off_est = lags(idx)/params.tpe_interp;
    phase_est = rad2deg(angle(r(idx)));
    
    % Align and estimate SNR
    % NOTE: this probably won't be very useful because there are multiple
    % signals (CDMA). Could modify the SNR estimator so it uses the peak
    % and noise floor out of the correlator.
    tpe_samps = rx_samps_interp(lags(idx)+1:end);
    tpe_samps = resample(tpe_samps, 1, params.tpe_interp);
    if length(tpe_samps) < length(ref_samps)
        n_pad = length(ref_samps) - length(tpe_samps);
        tpe_samps = [tpe_samps; zeros(n_pad, 1)];
    end
    idxs = 1:length(ref_samps);
    snr_samps = tpe_samps(idxs)*exp(-1j*deg2rad(phase_est));
    snr_samps = snr_samps/std(snr_samps);
    noise_samps = snr_samps - ref_samps;
    snr_est = 10*log10(1/var(noise_samps));
    
    % NOTE: the phase error here is no longer what is expected due to the
    % follower's adjustments. We should also check the rest of the errors
    % vs. expected if debugging is required.
    exp_phase = leader.rx_phase_off - ...
                rad2deg(leader.time*leader.rx.act_freq*2*pi) + ...
                l2f_chans(i).phase + ...
                followers(i).tx_phase_off + ...
                rad2deg(followers(i).time*followers(i).tx.act_freq*2*pi);
    exp_phase = mod(exp_phase, 360);
    phase_est_err = phase_est - exp_phase;
    
    % Status
    fprintf('Frequency offset : N/A\n');
    fprintf('Frequency error  : N/A\n');
    fprintf('Time offset      : %.4f samples\n', n_off_est);
    fprintf('Time error       : N/A\n');
    fprintf('Phase offset     : %.4f degrees\n', phase_est);
    fprintf('Phase error      : %.4f degrees\n', phase_est_err);
    fprintf('SNR              : %.4f dB\n', snr_est);
    fprintf('SNR error        : N/A\n');
    
    % Save data
    leader.data.f2l.followers(i).n_off_est = n_off_est;
    leader.data.f2l.followers(i).phase_est = phase_est;
    leader.data.f2l.followers(i).phase_est_err = phase_est_err;
    leader.data.f2l.followers(i).snr_est = snr_est;

end

%% All nodes delay and update clocks

% Delay includes the waveform duration as well as any protocol delay or
% processing times
n_delay = params.gc_len + params.proc_delay_len;
t_delay = n_delay/params.fs;

% Status
fprintf('\nProcessing delay, all nodes update clocks by %g seconds\n', ...
        t_delay);

% Update time for all nodes. Note that the follower clocks already include
% their offsets from the initial time of arrival estimates.
leader.time = leader.time + t_delay;
for i = 1:n_followers
    followers(i).time = followers(i).time + t_delay;
end

%% Leader sends channel estimates to followers

% Here we assume there is an error-free link and send a single message
% containing all of the follower channels. The followers receive and decode
% these for use during the beamforming step.

% NOTE: another idea may to be send a Gold code followed by payload bits,
% so that the followers can then refine their channel estimates.

for i = 1:n_followers
    
    % Status
    if i == 1
        fprintf('\n');
    end
    fprintf('Leader sends channel estimate to follower %i\n', i);
    
    % Assume data is sent without error
    followers(i).data.f2l.phase_est = leader.data.f2l.followers(i).phase_est;
    
end

%% All nodes delay and update clocks

% Delay includes the waveform duration as well as any protocol delay or
% processing times
n_delay = params.gc_len + params.proc_delay_len;
t_delay = n_delay/params.fs;

% Status
fprintf('\nProcessing delay, all nodes update clocks by %g seconds\n', ...
        t_delay);

% Update time for all nodes. Note that the follower clocks already include
% their offsets from the initial time of arrival estimates.
leader.time = leader.time + t_delay;
for i = 1:n_followers
    followers(i).time = followers(i).time + t_delay;
end

%% Delay waiting for target to transmit

% NOTE: to make this more realistic we could wait for some duration to
% detect the target signal. If it's not detected, then return to the start
% (leader transmits to the followers).

% Status
fprintf(['\nWaiting for target to transmit, all nodes update clocks ', ...
         'by %g seconds\n'], params.target_t_delay);

% Update time for all nodes. Note that the follower's clock offset from the
% leader is removed here, since receiving from the target does not depend
% on that (but it must be added back in before we transmit again).
leader.time = leader.time + params.target_t_delay;
for i = 1:n_followers
    n_off_est = followers(i).data.l2f.n_off_est;
    followers(i).time = followers(i).time - n_off_est/params.fs + ...
        params.target_t_delay;
end

%% Target transmits to all nodes

% Status
fprintf('\nTarget transmits to all nodes\n');

% Generate and modulate bits
% NOTE: using a Gold code for now. May want to consider other sequences,
% include different modulations, etc. in the future.
gc_gen = gold_code_gen(params.leader.id, params.gc_len);
bb_mod = bpsk_mod(params.alpha, params.span, params.sps);
bits = gc_gen.gen_bits();
ref_samps = bb_mod.sim(bits);

% TX simulation
% The signal starts at the target's current time
tx_samps = target.tx.sim(ref_samps, target.time);

% Save data
target.data.tx_bits = bits;
target.data.ref_samps = ref_samps;
target.data.tx_samps = tx_samps;

%% Leader receives from target

% Status
fprintf('\nLeader receives from target\n');

% Input data
tx_samps = target.data.tx_samps;

% Channel simulation
chan_samps = t2l_chan.sim(tx_samps);

% RX simulation
% NOTE: signal starts at the leader's current time, since the output of the
% channel includes the fractional sample delay (so the receiver phase will
% be set correctly by the time the signal arrives).
rx_samps = leader.rx.sim(chan_samps, leader.time);

% Save data
leader.data.t2l.chan_samps = chan_samps;
leader.data.t2l.rx_samps = rx_samps;

%% Followers receive from target

% Input data
tx_samps = target.data.tx_samps;

for i = 1:n_followers
    
    % Status
    if i == 1
        fprintf('\n');
    end
    fprintf('Follower %i receives from target\n', i);
    
    % Channel simulation
    chan_samps = t2f_chans(i).sim(tx_samps);
    
    % RX simulation
    % NOTE: signal starts at the follower's current time, since the output
    % of the channel includes the fractional sample delay (so the receiver
    % phase will be set correctly by the time the signal arrives).
    rx_samps = followers(i).rx.sim(chan_samps, followers(i).time);
    
    % Save data
    followers(i).data.t2f.chan_samps = chan_samps;
    followers(i).data.t2f.rx_samps = rx_samps;
    
end

%% Leader estimates channel from target

% Status
fprintf('\nLeader estimates channel from target\n');

% Input data
rx_samps = leader.data.t2l.rx_samps;

% Carrier frequency offset estimation and correction
cfc = comm.CoarseFrequencyCompensator('Modulation', 'BPSK', ...
                                      'FrequencyResolution', params.cfc_freq_res, ...
                                      'SampleRate', params.fs);
[cfc_samps, cfc_est] = cfc(rx_samps);

% Correlate to estimate time and phase offsets
% NOTE: could tweak this to assume we only know part of the signal
ref_samps = target.data.ref_samps;
cfc_samps_interp = resample(cfc_samps, params.tpe_interp, 1);
ref_samps_interp = resample(ref_samps, params.tpe_interp, 1);
[r, lags] = xcorr(cfc_samps_interp, ref_samps_interp);
[v, idx] = max(abs(r).^2);
n_off_est = lags(idx)/params.tpe_interp;
phase_est = rad2deg(angle(r(idx)));

% Align and estimate SNR
tpe_samps = cfc_samps_interp(lags(idx)+1:end);
tpe_samps = resample(tpe_samps, 1, params.tpe_interp);
idxs = 1:length(ref_samps);
snr_samps = tpe_samps(idxs)*exp(-1j*deg2rad(phase_est));
snr_samps = snr_samps/std(snr_samps);
noise_samps = snr_samps - ref_samps;
snr_est = 10*log10(1/var(noise_samps));

% Compute errors vs. expected
cfc_est_err = cfc_est - (target.freq_off - leader.freq_off);
n_off_est_err = n_off_est - t2l_chan.n_delay;
exp_phase = leader.rx_phase_off - ...
            rad2deg(leader.time*leader.rx.act_freq*2*pi) + ...
            t2l_chan.phase + ...
            target.tx_phase_off + ...
            rad2deg(target.time*target.tx.act_freq*2*pi);
exp_phase = mod(exp_phase, 360);
phase_est_err = phase_est - exp_phase;
snr_est_err = snr_est - (target.pow - t2l_chan.path_loss - leader.rx.noise_pow);

% Status
fprintf('Frequency offset : %.4f Hz\n', cfc_est)
fprintf('Frequency error  : %.4f Hz\n', cfc_est_err);
fprintf('Time offset      : %.4f samples\n', n_off_est);
fprintf('Time error       : %.4f samples\n', n_off_est_err);
fprintf('Phase offset     : %.4f degrees\n', phase_est);
fprintf('Phase error      : %.4f degrees\n', phase_est_err);
fprintf('SNR              : %.4f dB\n', snr_est);
fprintf('SNR error        : %.4f dB\n', snr_est_err);

% Save data
leader.data.t2l.cfc_est = cfc_est;
leader.data.t2l.cfc_est_err = cfc_est_err;
leader.data.t2l.n_off_est = n_off_est;
leader.data.t2l.n_off_est_err = n_off_est_err;
leader.data.t2l.phase_est = phase_est;
leader.data.t2l.phase_est_err = phase_est_err;
leader.data.t2l.snr_est = snr_est;
leader.data.t2l.snr_est_err = snr_est_err;

%% Followers estimate channels from target

for i = 1:n_followers
    
    % Status
    fprintf('\nFollower %i estimates channel from target\n', i);
    
    % Input data
    rx_samps = followers(i).data.t2f.rx_samps;
    
    % Carrier frequency offset estimation and correction
    cfc = comm.CoarseFrequencyCompensator('Modulation', 'BPSK', ...
                                          'FrequencyResolution', params.cfc_freq_res, ...
                                          'SampleRate', params.fs);
    [cfc_samps, cfc_est] = cfc(rx_samps);
    
    % Correlate to estimate time and phase offsets
    % NOTE: could tweak this to assume we only know part of the signal
    ref_samps = target.data.ref_samps;
    cfc_samps_interp = resample(cfc_samps, params.tpe_interp, 1);
    ref_samps_interp = resample(ref_samps, params.tpe_interp, 1);
    [r, lags] = xcorr(cfc_samps_interp, ref_samps_interp);
    [v, idx] = max(abs(r).^2);
    n_off_est = lags(idx)/params.tpe_interp;
    phase_est = rad2deg(angle(r(idx)));
    
    % Adjust phase estimate based on the protocol delay (since the time we
    % estimated the frequency offset from the leader) and the delay of the
    % target signal arriving
    f_corr = followers(i).data.l2f.cfc_est;
    t_err = (n_off_est - t2f_chans(i).n_delay)/params.fs;
    t_total = 6e-3 + (params.target_t_delay + t_err);
    phase_est = phase_est - rad2deg(2*pi*f_corr*t_total);

    % Align and estimate SNR
    tpe_samps = cfc_samps_interp(lags(idx)+1:end);
    tpe_samps = resample(tpe_samps, 1, params.tpe_interp);
    idxs = 1:length(ref_samps);
    snr_samps = tpe_samps(idxs)*exp(-1j*deg2rad(phase_est));
    snr_samps = snr_samps/std(snr_samps);
    noise_samps = snr_samps - ref_samps;
    snr_est = 10*log10(1/var(noise_samps));
    
    % Compute errors vs. expected
    % NOTE: the phase error here is no longer what is expected due to the
    % follower's adjustments.
    cfc_est_err = cfc_est - (target.freq_off - followers(i).freq_off);
    n_off_est_err = n_off_est - t2f_chans(i).n_delay;
    exp_phase = target.tx_phase_off + ...
                rad2deg(target.time*target.tx.act_freq*2*pi) + ...
                t2f_chans(i).phase + ...
                followers(i).rx_phase_off - ...
                rad2deg(followers(i).time*followers(i).rx.act_freq*2*pi);
    exp_phase = mod(exp_phase, 360);
    phase_est_err = phase_est - exp_phase;
    snr_est_err = snr_est - (target.pow - t2f_chans(i).path_loss - followers(i).rx.noise_pow);
    
    % Status
    fprintf('Frequency offset : %.4f Hz\n', cfc_est)
    fprintf('Frequency error  : %.4f Hz\n', cfc_est_err);
    fprintf('Time offset      : %.4f samples\n', n_off_est);
    fprintf('Time error       : %.4f samples\n', n_off_est_err);
    fprintf('Phase offset     : %.4f degrees\n', phase_est);
    fprintf('Phase error      : %.4f degrees\n', phase_est_err);
    fprintf('SNR              : %.4f dB\n', snr_est);
    fprintf('SNR error        : %.4f dB\n', snr_est_err);
    
    % Save data
    followers(i).data.t2f.cfc_est = cfc_est;
    followers(i).data.t2f.cfc_est_err = cfc_est_err;
    followers(i).data.t2f.n_off_est = n_off_est;
    followers(i).data.t2f.n_off_est_err = n_off_est_err;
    followers(i).data.t2f.phase_est = phase_est;
    followers(i).data.t2f.phase_est_err = phase_est_err;
    followers(i).data.t2f.snr_est = snr_est;
    followers(i).data.t2f.snr_est_err = snr_est_err;
    
end

%% All nodes delay and update clocks

% Delay includes the waveform duration as well as any protocol delay or
% processing times
% NOTE: if the target signal is changed, then update the waveform duration
n_delay = params.gc_len + params.proc_delay_len;
t_delay = n_delay/params.fs;

% Status
fprintf('\nProcessing delay, all nodes update clocks by %g seconds\n', ...
        t_delay);

% Update time for all nodes. Note that we add back in the follower's clock
% offset since it is transmitting based on the estimated time of arrival
% from the initial leader's signal. Also note that the target's time is
% updated here, but that doesn't really matter, and we can add a random
% offset to demonstrate this.
% NOTE: basing the next transmit time off of the target signal's estimated
% time of arrival will break things. Because of this we remove the delay
% waiting for the target signal before updating, so we end up at the next
% slot time relative to the leader's clock.
leader.time = leader.time - params.target_t_delay + t_delay;
for i = 1:n_followers
    n_off_est = followers(i).data.l2f.n_off_est;
    followers(i).time = followers(i).time - params.target_t_delay + ...
        n_off_est/params.fs + t_delay;
end
target.time = target.time + randn*10e-9 + t_delay;

%% All nodes transmit beamforming signal

% Generate and modulate bits
% NOTE: using the leader's Gold code for now. May want to consider other
% sequences, include different modulations, etc. in the future.
gc_gen = gold_code_gen(params.leader.id, params.gc_len);
bb_mod = bpsk_mod(params.alpha, params.span, params.sps);
bits = gc_gen.gen_bits();
ref_samps = bb_mod.sim(bits);

% Status
fprintf('\nLeader transmits beamforming signal to target\n');

% Phase shift for beamforming
phase_shift = -leader.data.t2l.phase_est;
ref_samps_phase = ref_samps*exp(1j*deg2rad(phase_shift));

% TX simulation
% The signal starts at the leader's current time
tx_samps = leader.tx.sim(ref_samps_phase, leader.time);

% Save data
leader.data.l2t.bits = bits;
leader.data.l2t.ref_samps = ref_samps;
leader.data.l2t.tx_samps = tx_samps;

for i = 1:n_followers
    
    % Status
    fprintf('Follower %i transmits beamforming signal to target\n', i);
        
    % Apply estimated frequency offset
    cfc_est = followers(i).data.l2f.cfc_est;
    n = 0:(length(ref_samps) - 1);
    ref_samps_freq = ref_samps.*exp(1j*2*pi*(cfc_est/params.fs)*n');

    % Apply phase offset based on the protocol delay (since the time we
    % estimated the frequency offset from the leader)
    ref_samps_freq = ref_samps_freq*exp(1j*2*pi*cfc_est*8e-3);

    % NOTE: we could adjust for all of the phase terms due to frequency
    % offset here, it would be by 12 ms + target delay. The 12 ms is
    % because of the signs below (0, -2, +6, +8).

    % Phase shift for beamforming
    phase_shift = followers(i).data.l2f.phase_est - ...
                  followers(i).data.f2l.phase_est - ...
                  followers(i).data.t2f.phase_est;
    ref_samps_phase = ref_samps_freq*exp(1j*deg2rad(phase_shift));
    
    % TX simulation
    % The signal starts at the follower's current time
    tx_samps = followers(i).tx.sim(ref_samps_phase, followers(i).time);

    % Save data
    followers(i).data.f2t.bits = bits;
    followers(i).data.f2t.ref_samps = ref_samps;
    followers(i).data.f2t.tx_samps = tx_samps;
    
end

%% Target receives from all nodes

% Status
fprintf('\nTarget receives from leader\n');

% Input data
tx_samps = leader.data.l2t.tx_samps;

% Channel simulation
chan_samps = t2l_chan.sim(tx_samps);

for i = 1:n_followers
    
    % Status
    fprintf('Target receives from follower %i\n', i);
    
    % Input data
    tx_samps = followers(i).data.f2t.tx_samps;
    
    % Channel simulation (save each signal independently)
    chan_samps = [chan_samps, t2f_chans(i).sim(tx_samps)];
    
end

% Target will receive summation of all signals
sum_chan_samps = sum(chan_samps, 2);

% RX simulation
% NOTE: signal starts at the target's current time, since the output of the
% channel includes the fractional sample delay (so the receiver phase will
% be set correctly by the time the signal arrives).
rx_samps = target.rx.sim(sum_chan_samps, target.time);

% Save data
target.data.chan_samps = chan_samps;
target.data.rx_samps = rx_samps;

%% Target receiver processing

% Status
fprintf('\nTarget receiver processing\n');

% NOTE: this is only with beamforming. For debugging we may want to compare
% the result to each individual signal, which would need to be run through
% the receiver to add noise and etc.

% Input data
rx_samps = target.data.rx_samps;

% Carrier frequency offset estimation and correction
cfc = comm.CoarseFrequencyCompensator('Modulation', 'BPSK', ...
                                      'FrequencyResolution', params.cfc_freq_res, ...
                                      'SampleRate', params.fs);
[cfc_samps, cfc_est] = cfc(rx_samps);

% Correlate to estimate time and phase offsets
% NOTE: could try assuming that we only know part of the signal
ref_samps = leader.data.l2t.ref_samps;
cfc_samps_interp = resample(cfc_samps, params.tpe_interp, 1);
ref_samps_interp = resample(ref_samps, params.tpe_interp, 1);
[r, lags] = xcorr(cfc_samps_interp, ref_samps_interp);
[v, idx] = max(abs(r).^2);
n_off_est = lags(idx)/params.tpe_interp;
phase_est = rad2deg(angle(r(idx)));

% Align and estimate SNR
tpe_samps = cfc_samps_interp(lags(idx)+1:end);
tpe_samps = resample(tpe_samps, 1, params.tpe_interp);
idxs = 1:length(ref_samps);
snr_samps = tpe_samps(idxs)*exp(-1j*deg2rad(phase_est));
snr_samps = snr_samps/std(snr_samps);
noise_samps = snr_samps - ref_samps;
snr_est = 10*log10(1/var(noise_samps));

% Matched filter
rrc_filt = comm.RaisedCosineReceiveFilter(...
                'RolloffFactor', params.alpha, ...
                'FilterSpanInSymbols', params.span, ...
                'InputSamplesPerSymbol', params.sps, ...
                'DecimationFactor', params.sps, ...
                'DecimationOffset', 0);
rx_syms = rrc_filt(snr_samps);

% Remove the delay due to both TX and RX filters
n_delay = 1 + 2*(rrc_filt.order/2)/8;
rx_syms = rx_syms(n_delay+1:end);

% Make bit decisions and compute BER
rx_bits = rx_syms < 0;
n_errs = sum(xor(leader.data.l2t.bits, rx_bits));
ber = n_errs/length(rx_bits);

% NOTE: could compute errors vs. expected if debugging is required

% Status
fprintf('Frequency offset : %.4f Hz\n', cfc_est)
fprintf('Time offset      : %.4f samples\n', n_off_est);
fprintf('Phase offset     : %.4f degrees\n', phase_est);
fprintf('SNR              : %.4f dB\n', snr_est);
fprintf('Bit errors       : %i\n', n_errs);
fprintf('BER              : %.4f\n', ber);

% Save data
target.data.cfc_est = cfc_est;
target.data.n_off_est = n_off_est;
target.data.phase_est = phase_est;
target.data.snr_est = snr_est;
target.data.rx_syms = rx_syms;
target.data.rx_bits = rx_bits;
target.data.n_errs = n_errs;
target.data.ber = ber;

%{
% Plot the received constellation
figure
plot(rx_syms, 'o')
grid on
axis equal
%}

%% Evaluation of beamforming gain

% Ideal beamforming gain
% This is N^2, which assumes equal transmit powers and distances
ideal_bf_gain = 20*log10(1 + n_followers);

% Theoretical maximum received power (dBm)
% This is based on the actual transmit powers and path losses
theory_max_bf_power = sqrt(10^(((leader.pow - t2l_chan.path_loss) - 30)/10));
for i = 1:n_followers
    theory_max_bf_power = theory_max_bf_power + ...
        sqrt(10^(((followers(i).pow - t2f_chans(i).path_loss) - 30)/10));
end
theory_max_bf_power = 20*log10(theory_max_bf_power) + 30;

% Theoretical maximum beamforming gains (dB)
% Average and relative to each node
theory_bf_gains = zeros(1 + n_followers, 1);
theory_bf_gains(1) = theory_max_bf_power - (leader.pow - t2l_chan.path_loss);
for i = 2:(n_followers + 1)
    theory_bf_gains(i) = theory_max_bf_power - (followers(i-1).pow - t2f_chans(i-1).path_loss);
end
avg_theory_bf_gain = 10*log10(mean(10.^(theory_bf_gains/10)));

% Actual received power (dBm)
% NOTE: this includes (signal + noise) powers, so will only be accurate at
% high SNRs. Should also look at the SNR gain to verify.
act_bf_power = 20*log10(sqrt(var(sum_chan_samps))) + 30;

% Actual beamforming gains (dB)
% Relative to each node and the average
% NOTE: this includes (signal + noise) powers, so will only be accurate at
% high SNRs. Should also look at the SNR gain to verify.
act_bf_gains = zeros(1 + n_followers, 1);
act_bf_gains(1) = act_bf_power - (20*log10(std(chan_samps(:, 1))) + 30);
for i = 2:(n_followers + 1)
    act_bf_gains(i) = act_bf_power - (20*log10(std(chan_samps(:, i))) + 30);
end
avg_act_bf_gain = 10*log10(mean(10.^(act_bf_gains/10)));

% Status
fprintf('Ideal beamforming gain (N^2) : %.4f dB\n', ideal_bf_gain);
fprintf('Theoretical max. RX power    : %.4f dBm\n', theory_max_bf_power);
fprintf('Theoretical max. beamforming gains\n');
fprintf('Average                      : %.4f dB\n', avg_theory_bf_gain);
fprintf('Relative to leader           : %.4f dB\n', theory_bf_gains(1));
for i = 2:(n_followers + 1)
    fprintf('Relative to follower %i       : %.4f dB\n', (i-1), theory_bf_gains(i));
end
fprintf('Actual RX power              : %.4f dBm\n', act_bf_power);
fprintf('Actual beamforming gains\n');
fprintf('Average                      : %.4f dB\n', 10*log10(mean(10.^(act_bf_gains/10))));
fprintf('Relative to leader           : %.4f dB\n', act_bf_gains(1));
for i = 2:(n_followers + 1)
    fprintf('Relative to follower %i       : %.4f dB\n', (i-1), act_bf_gains(i));
end
%% Plots
%{
% Individual and combined signals
figure
for i = 1:(n_followers + 1)
        
    subplot(2,1,1)
    hold on
    plot(real(target.data.chan_samps(:, i)))
    hold off
    grid on
    ylabel('I')
    
    subplot(2,1,2)
    hold on
    plot(imag(target.data.chan_samps(:, i)))
    hold off
    grid on
    ylabel('Q')

end
subplot(2,1,1)
hold on
plot(real(sum_chan_samps), 'k')
hold off
subplot(2,1,2)
hold on
plot(imag(sum_chan_samps), 'k')
hold off

% Individual constellations
figure
hold on
for i = 1:(n_followers + 1)
    plot(target.data.chan_samps(:, i), 'o')
end
hold off
grid on
axis equal

% Individual constellations with frequency offsets removed

r = target.data.chan_samps(:, 1);
r = r/max(abs(r));
s = exp(-1j*2*pi*(leader.freq_off/params.fs)*(0:length(r)-1)');

figure
hold on
plot(r.*s, '.')

r = target.data.chan_samps(:, 2);
r = r/max(abs(r));
s = exp(-1j*2*pi*(leader.freq_off/params.fs)*(0:length(r)-1)');

plot(r.*s, '.')

r = target.data.chan_samps(:, 3);
r = r/max(abs(r));
s = exp(-1j*2*pi*(leader.freq_off/params.fs)*(0:length(r)-1)');

plot(r.*s, '.')

hold off
grid on
xlim([-1,1])
ylim([-1,1])
%}