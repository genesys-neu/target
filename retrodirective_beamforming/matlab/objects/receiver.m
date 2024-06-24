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

classdef receiver < handle
    
    %RECEIVER RF receiver model
    %   Models noise, frequency offset, and phase.
    
    properties
        fs
        noise_fig
        freq
        freq_off
        act_freq
        phase_off
        noise_pow
        noise
        sig_pow
        snr
    end
    
    methods
        
        function obj = receiver(fs, noise_fig, freq, freq_off, phase_off)
            
            % Copy inputs
            %   fs        - sample rate (SPS)
            %   noise_fig - noise figure (dB)
            %   freq      - carrier frequency (Hz)
            %   freq_off  - frequency offset (Hz)
            %   phase_off - phase offset (degrees)
            obj.fs =  fs;
            obj.freq = freq;
            obj.noise_fig = noise_fig;
            obj.freq_off = freq_off;
            obj.phase_off = phase_off;
            
            % Noise power (dBm)
            k = physconst('boltzmann');
            t = 290;
            obj.noise_pow = 10*log10(k*t*obj.fs) + 30 + obj.noise_fig;
            
            % Actual frequency (Hz)
            obj.act_freq = obj.freq + obj.freq_off;
            
        end
        
        function y = sim(obj, x, t)
            
            % Inputs
            %   x - signal (column vector)
            %   t - current time (seconds)
                        
            % Make sure input is a column vector
            x = x(:);
            
            % Add noise
            n = length(x);
            noise_std = sqrt(10^((obj.noise_pow - 30)/10)/2);
            obj.noise = noise_std*(randn(n, 1) + 1j*randn(n, 1));
            x_noise = x + obj.noise;
            
            % Save signal power (dBm)
            % NOTE: signal power is computed using the entire input,
            % including any zero padding.
            obj.sig_pow = 10*log10(var(x)) + 30;
            
            % Save SNR (dB)
            obj.snr = (obj.sig_pow - 30) - 10*log10(var(obj.noise));

            % Phase offset (initial offset + offset at current time)
            curr_phase = obj.phase_off - rad2deg(t*obj.act_freq*2*pi);
            x_phase = x_noise*exp(1j*deg2rad(curr_phase));
            
            % Frequency offset
            n = 0:(length(x_phase) - 1);
            x_freq = x_phase.*exp(-1j*2*pi*(obj.freq_off/obj.fs)*n');
            
            % Output
            y = x_freq;
            
        end
        
    end
    
end