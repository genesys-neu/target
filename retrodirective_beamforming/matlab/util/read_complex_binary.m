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

function x = read_complex_binary(filename, n, offset)
%READ_COMPLEX_BINARY Read samples from complex float (fc32) binary file

%   x        : signal (floating point)
%   filename : input file
%   n        : number of samples to read (default: all)
%   offset   : number of samples to skip before reading (default: 0)

% Set default arguments
if nargin < 2 
    n = Inf;
    offset = 0;
elseif nargin < 3
    offset = 0;
end

% Open the file
fid = fopen(filename, 'rb');
if fid < 0
    warning('Could not open file');
end

% Skip samples if necessary (each sample is 8 bytes)
if offset > 0
    status = fseek(fid, offset*8, 'bof');
    if status < 0
        warning('Could not skip requested samples from beginning of file');
    end
end

% Read from the file
[y, count] = fread(fid, 2*n, 'float');
if count < 2*n && n ~= Inf
    warning('Could not read all requested samples from file')
end

% Close the file
fid = fclose(fid);
if fid < 0
    warning('Could not close file');
end

% Convert to complex floating point
x = y(1:2:end) + 1j*y(2:2:end);

end