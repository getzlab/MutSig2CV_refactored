function ni = count_overlaps_fast2(keys, pat_idx)
% COUNT_OVERLAPS_FAST2 - MATLAB implementation of fast overlap counting
% This replaces the MEX-based count_overlaps_fast2 function for compatibility
%
% ni = count_overlaps_fast2(keys, pat_idx)
% Returns a patient-patient matrix where ni(i,j) = number of overlapping mutations
% between patient i and patient j
%
% Input:
%   keys - sorted mutation keys (e.g., chr*1e9+start)
%   pat_idx - patient indices corresponding to each mutation
%
% Output:
%   ni - npat x npat matrix of overlap counts (diagonal = 0)

% Input validation
if nargin ~= 2
    error('count_overlaps_fast2 requires exactly 2 arguments');
end

% Ensure inputs are column vectors
keys = keys(:);
pat_idx = pat_idx(:);

% Check that inputs have same length
if length(keys) ~= length(pat_idx)
    error('keys and pat_idx must have the same length');
end

% Remove any NaN or Inf values
valid = ~isnan(keys) & ~isnan(pat_idx) & ~isinf(keys) & ~isinf(pat_idx);
keys = keys(valid);
pat_idx = pat_idx(valid);

% Check that keys are sorted
if any(diff(keys) < 0)
    error('keys must be sorted in ascending order');
end

% Check that pat_idx values are valid
if any(pat_idx < 1)
    error('pat_idx must be 1 or greater');
end

% Find number of patients
npat = max(pat_idx);

% Initialize output matrix
ni = zeros(npat, npat);

% Count overlaps
i = 1;
while i <= length(keys)
    j = i;
    keyi = keys(i);
    
    % Find all mutations with the same key
    while j + 1 <= length(keys) && keys(j + 1) == keyi
        j = j + 1;
    end
    
    % If there are multiple mutations with the same key, count overlaps
    if j > i
        for a = i:j
            pa = pat_idx(a);
            for b = i:j
                pb = pat_idx(b);
                ni(pb, pa) = ni(pb, pa) + 1;
            end
        end
    end
    
    i = j + 1;
end

% Clear the diagonal (no self-overlaps)
for p = 1:npat
    ni(p, p) = 0;
end 