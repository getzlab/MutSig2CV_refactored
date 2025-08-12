function H = hist2d_fast(x, y, xmin, xmax, ymin, ymax)
% HIST2D_FAST - MATLAB implementation of 2D histogram
% This replaces the MEX-based hist2d_fast function for compatibility
%
% H = hist2d_fast(x, y, xmin, xmax, ymin, ymax)
% Creates a 2D histogram of data points (x,y) with specified ranges

% Input validation
if nargin ~= 6
    error('hist2d_fast requires exactly 6 arguments');
end

% Ensure inputs are column vectors
x = x(:);
y = y(:);

% Remove any NaN or Inf values
valid = ~isnan(x) & ~isnan(y) & ~isinf(x) & ~isinf(y);
x = x(valid);
y = y(valid);

% Calculate bin edges
nx = xmax - xmin + 1;
ny = ymax - ymin + 1;

% Initialize output matrix
H = zeros(nx, ny);

% Create bin indices
x_bins = x - xmin + 1;
y_bins = y - ymin + 1;

% Filter out-of-bounds values
in_bounds = x_bins >= 1 & x_bins <= nx & y_bins >= 1 & y_bins <= ny;
x_bins = x_bins(in_bounds);
y_bins = y_bins(in_bounds);

% Use accumarray for fast counting
if ~isempty(x_bins)
    H = accumarray([x_bins, y_bins], 1, [nx, ny]);
end 