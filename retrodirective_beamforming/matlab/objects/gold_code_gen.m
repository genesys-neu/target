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

classdef gold_code_gen
    
    %GOLD_CODE_GEN Gold code generator
    %   Generates samples of a gold code for channel sounding.
    
    properties
        id
        n
        gc
    end
    
    methods
        
        function obj = gold_code_gen(id, n)
                                    
            % Check inputs
            if (id < -2) || (id > (2^10 - 2))
                error('ID is out of range');
            end
            if (n < 1) || (n > (2^10 - 1))
                error('n is out of range');
            end
            
            % Copy inputs
            %   id - sequence index
            %   n  - samples per frame
            obj.id = id;
            obj.n = n;
                        
            % Make the gold code generator
            % This uses the preferred pair with polynomial degree 10 and a
            % sequence length of 1023.
            obj.gc = comm.GoldSequence(...
                        'FirstPolynomial', [10 3 0], ...
                        'FirstInitialConditions', [0 0 0 0 0 0 0 0 0 1], ...
                        'SecondPolynomial', [10 8 3 2 0], ...
                        'SecondInitialConditions', [0 0 0 0 0 0 0 0 0 1], ...
                        'Index', id,  ...
                        'SamplesPerFrame', n);
            
        end
        
        function reset(obj)
            
            obj.gc.reset();
            
        end
        
        function bits = gen_bits(obj)
            
            % Reset so we always generate the same pattern
            obj.reset();
            
            % Generate bits
            bits = obj.gc();
            
        end
        
    end
    
end