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

classdef node < handle
    
    %NODE Node that has a transmitter and receiver
    %   Models a node that has coordinates and an RF transmitter/receiver.
    
    properties
        id
        xyz
        fs
        time
        freq
        freq_off
        tx_phase_off
        rx_phase_off
        noise_fig
        pow
        tx
        rx
        data
    end
    
    methods
        
        function obj = node(params)
            
            % Unpack parameters
            %   id           - identifier
            %   xyz          - vector of cartesian coordinates (meters)
            %   fs           - sample rate (SPS)
            %   time         - clock time (seconds)
            %   freq         - carrier frequency (Hz)
            %   freq_off     - frequency offset (Hz)
            %   tx_phase_off - transmitter phase offset (degrees)
            %   rx_phase_off - receiver phase offset (degrees)
            %   noise_fig    - noise figure (dB)
            %   pow          - output power (dBm)
            obj.id = params.id;
            obj.xyz = params.xyz(:)';
            obj.fs = params.fs;
            obj.time = params.time;
            obj.freq = params.freq;
            obj.freq_off = params.freq_off;
            obj.tx_phase_off = params.tx_phase_off;
            obj.rx_phase_off = params.rx_phase_off;
            obj.noise_fig = params.noise_fig;
            obj.pow = params.pow;
            
            % Note that the transmitter and receiver share sample rate and
            % frequency error, but have independent phases.
            
            % Initialize transmitter
            obj.tx = transmitter(obj.fs, ...
                                 obj.freq, ...
                                 obj.freq_off, ...
                                 obj.tx_phase_off, ...
                                 obj.pow);
                             
            % Initialize receiver
            obj.rx = receiver(obj.fs, ...
                              obj.noise_fig, ...
                              obj.freq, ...
                              obj.freq_off, ...
                              obj.rx_phase_off);
                          
            % Make a struct for saving this node's data
            obj.data = struct;
            
        end
        
    end
    
end