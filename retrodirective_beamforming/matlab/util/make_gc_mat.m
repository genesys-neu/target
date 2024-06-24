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

% Makes Gold codes and saves to file for experiments with SDRs

clear

% Make codes for up to 10 nodes, 1000 bits each
node_ids = 0:9;
n_bits = 1000;

n_nodes = length(node_ids);
b_mat = zeros(n_nodes, n_bits);
for i = 1:n_nodes
    
    gc = gold_code_gen(node_ids(i), n_bits);
    b_mat(i, :) = gc.gen_bits();
    
end

% Save to .mat file
save('gc.mat', 'b_mat');