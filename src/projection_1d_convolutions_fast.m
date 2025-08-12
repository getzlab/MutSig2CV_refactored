function [pmax, pmin] = projection_1d_convolutions_fast(Sdeg, Pdeg, score_obs, numbins, H, newH)
% PROJECTION_1D_CONVOLUTIONS_FAST - MATLAB implementation of 1D convolution projection
% This replaces the MEX-based projection_1d_convolutions_fast function for compatibility
%
% [pmax, pmin] = projection_1d_convolutions_fast(Sdeg, Pdeg, score_obs, numbins, H, newH)
%
% Input:
%   Sdeg - score of each 1D degree (np,ncat+1): MUST BE INTEGERS
%   Pdeg - probability of each 1D degree (np,ncat+1)
%   score_obs - observed score to calculate p-value for
%   numbins - typically score_obs + ~10
%   H - externally allocated matrix (numbins,1)
%   newH - externally allocated matrix (numbins,ncols)
%
% Returns:
%   pmax - upper bound on p-value, i.e. the standard p-value usually calculated
%   pmin - lower bound on p-value, going "one discrete step further"

% Input validation
if nargin ~= 6
    error('projection_1d_convolutions_fast requires exactly 6 arguments');
end

% Get dimensions
[np, ncat_plus_1] = size(Sdeg);
ncat = ncat_plus_1 - 1;
ncols = ncat + 1;

% Validate input sizes
if size(Pdeg, 1) ~= np || size(Pdeg, 2) ~= ncat_plus_1
    error('Pdeg should be same size as Sdeg');
end

if size(H, 1) ~= numbins || size(H, 2) ~= 1
    error('H size should be (numbins,1)');
end

if size(newH, 1) ~= numbins || size(newH, 2) ~= ncols
    error('newH size should be (numbins,ncols)');
end

% Initialize H: all probability in first bin
H(:) = 0;
H(1) = 1;

% Sequential convolution
for p = 1:np
    newH(:) = 0;
    
    for d = 1:ncat_plus_1
        odouble = Sdeg(p, d);
        o = round(odouble); % round to nearest integer
        
        if o < numbins
            idx2 = (d-1) * numbins + o + 1;
            for k = 1:(numbins - o)
                newH(idx2 + k - 1) = Pdeg(p, d) * H(k);
            end
        end
    end
    
    % Sum newH across columns to get H
    for i = 1:numbins
        H(i) = 0;
        for j = 1:ncols
            idx = i + (j-1) * numbins;
            H(i) = H(i) + newH(idx);
        end
    end
end

% Calculate p-value
pmax = 1;
pmin = 1;

for i = 1:numbins
    if H(i) > 0
        pmin = pmin - H(i);
        if i > score_obs
            break;
        end
        pmax = pmax - H(i);
    end
end

pmax = max(0, pmax);
pmin = max(0, pmin); 