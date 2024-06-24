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

classdef channel < handle
    
    %CHANNEL RF channel model
    %   Models time delay, phase shift, and path loss between two nodes.
    
    properties
        src
        dst
        fs
        freq
        lambda
        path_len
        delay
        phase
        path_loss
        n_delay
        delay_filt
    end
    
    methods
        
        function obj = channel(src, dst)
            
            % Copy inputs
            %   src - source node
            %   dst - destination node
            obj.src = src;
            obj.dst = dst;
            
            % Make sure sample rate and carrier frequency match
            if obj.src.fs ~= obj.dst.fs
                error('Sample rate mismatch');
            end
            if obj.src.freq ~= obj.dst.freq
                error('Carrier frequency mismatch');
            end
            obj.fs = obj.src.fs;
            obj.freq = obj.src.freq;
            
            % Path length (m)
            obj.path_len = pdist2(obj.src.xyz, obj.dst.xyz);
                                    
            % Speed of light (m/s)
            c = physconst('lightspeed');
            
            % Wavelength (m)
            obj.lambda = c/obj.freq;
            
            % Time delay (s)
            obj.delay = obj.path_len/c;
            
            % Phase shift (deg)
            lambdas = obj.path_len/obj.lambda;
            obj.phase = rad2deg(2*pi*(lambdas - floor(lambdas)));
            
            % Path loss (dB)
            obj.path_loss = 20*log10(4*pi*obj.path_len/obj.lambda);
            
            % Fractional delay filter
            obj.n_delay = obj.delay*obj.fs;
            frac_delay = obj.n_delay - floor(obj.n_delay);
            h = designFracDelayFIR(frac_delay);
            obj.delay_filt = dsp.FIRFilter(h);
            
        end
        
        function reset(obj)
            
            obj.delay_filt.reset();
            
        end
        
        function y = sim(obj, x)
            
            % NOTE: we could keep a buffer of samples so we can call this over
            % and over again without having to reset anything.
            obj.reset()
            
            % Make sure input is a column vector
            x = x(:);
            
            % Time delay on the baseband samples. First do the integer
            % part, then the fractional part. The input to the fractional
            % delay filter is padded so that its group delay can be removed
            % afterwards.
            int_delay = floor(obj.n_delay);
            x = [zeros(int_delay, 1); x];
            n_pad = floor(obj.delay_filt.order/2);
            x_pad = [x; zeros(n_pad, 1)];
            x_delay = obj.delay_filt(x_pad);
            x_delay(1:n_pad) = [];
            
            % Phase shift on the RF carrier. This is simulated since there
            % are many carrier cycles for each baseband sample (so this
            % wouldn't need to be done if we simulated the delay on the
            % carrier samples, but that's a really high sample rate).
            x_shift = x_delay*exp(1j*deg2rad(obj.phase));
            
            % Path loss
            x_loss = x_shift*10^(-obj.path_loss/20);
            
            % Output
            y = x_loss;
            
        end
        
    end
    
end