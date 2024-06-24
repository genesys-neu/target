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

classdef bpsk_mod < handle
    
    %BPSK_MOD BPSK modulator
    %   Takes input bits and generates BPSK modulated samples.
    
    properties
        alpha
        span
        sps
        rrc_filt
        n_pad
    end
    
    methods
        
        function obj = bpsk_mod(alpha, span, sps)

            % Copy inputs
            %   alpha - RRC roll-off factor
            %   span  - RRC span in symbols
            %   sps   - samples per symbol
            obj.alpha = alpha;
            obj.span = span;
            obj.sps = sps;
            
            % RRC filter
            obj.rrc_filt = comm.RaisedCosineTransmitFilter(...
                                'RolloffFactor', obj.alpha, ...
                                'FilterSpanInSymbols', obj.span, ...
                                'OutputSamplesPerSymbol', obj.sps);
            
            % Symbol padding to handle filter delay and transients
            obj.n_pad = 2*ceil(floor(obj.rrc_filt.order/2)/obj.sps);
            
        end
        
        function reset(obj)
           
            obj.rrc_filt.reset();
            
        end
        
        function x = sim(obj, bits)

            % NOTE: we could keep a buffer of samples so we can call this over
            % and over again without having to reset anything.
            obj.reset()
            
            % Make sure input is a column vector
            bits = bits(:);

            % Map bits to symbols
            syms = 1 - 2*bits;
            
            % Pad to account for filter delay and transients
            syms_pad = [0; syms; zeros(obj.n_pad, 1)];
            
            % Upsample and filter
            samps = obj.rrc_filt(syms_pad);
            
            % Normalize to unit power
            % NOTE: this also includes the padding
            samps_norm = samps/std(samps);
            
            % Output
            x = samps_norm;
                        
        end
        
    end
    
end