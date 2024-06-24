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

classdef transmitter < handle
    
    %TRANSMITTER RF transmitter model
    %   Models frequency offset, phase, and gain.
    
    properties
        fs
        freq
        freq_off
        act_freq
        phase_off
        pow
    end
    
    methods
        
        function obj = transmitter(fs, freq, freq_off, phase_off, pow)
            
            % Copy inputs
            %   fs        - sample rate (SPS)
            %   freq      - carrier frequency (Hz)
            %   freq_off  - frequency offset (Hz)
            %   phase_off - phase offset (degrees)
            %   pow       - output power (dBm)
            obj.fs =  fs;
            obj.freq = freq;
            obj.freq_off = freq_off;
            obj.phase_off = phase_off;
            obj.pow = pow;
            
            % Actual frequency (Hz)
            obj.act_freq = obj.freq + obj.freq_off;
            
        end
        
        function y = sim(obj, x, t)
            
            % Inputs
            %   x - signal (column vector)
            %   t - current time (seconds)
                        
            % Make sure input is a column vector
            x = x(:);
            
            % Phase offset (initial offset + offset at current time)
            curr_phase = obj.phase_off + rad2deg(t*obj.act_freq*2*pi);
            x_phase = x*exp(1j*deg2rad(curr_phase));
            
            % Frequency offset            
            n = 0:(length(x_phase) - 1);
            x_freq = x_phase.*exp(1j*2*pi*(obj.freq_off/obj.fs)*n');
            
            % Apply gain; output is in W
            var_out = 10^((obj.pow - 30)/10);
            x_gain = x_freq*sqrt(var_out)/std(x_freq);
            
            % Output
            y = x_gain;
            
        end
        
    end
    
end