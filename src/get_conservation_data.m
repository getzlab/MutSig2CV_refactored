function conservation_data = get_conservation_data(chr, start_pos, end_pos, conservation_fwb_file, FixedWidthBinary_jar_file)
%
% GET_CONSERVATION_DATA - Extract conservation data from FWB file
%
% USAGE:
%   conservation_data = get_conservation_data(chr, start_pos, end_pos, conservation_fwb_file, FixedWidthBinary_jar_file)
%
% INPUTS:
%   chr - chromosome identifier (can be string like 'chr1', '1', or numeric)
%   start_pos - start position (1-based)
%   end_pos - end position (1-based)
%   conservation_fwb_file - path to the conservation FWB file
%   FixedWidthBinary_jar_file - path to the FixedWidthBinary.jar file
%
% OUTPUTS:
%   conservation_data - structure containing:
%     .data - vector of conservation scores (NaN for missing data)
%     .positions - vector of genomic positions
%     .chr - chromosome number
%     .start - start position
%     .end - end position
%     .length - length of the region
%     .num_valid - number of valid (non-NaN) positions
%     .mean_conservation - mean conservation score (excluding NaN)
%     .min_conservation - minimum conservation score
%     .max_conservation - maximum conservation score
%
% EXAMPLE:
%   conservation_data = get_conservation_data('chr1', 1000000, 1000100, ...
%       'reference/conservation46.fwb', 'reference/FixedWidthBinary.jar');
%
% NOTES:
%   - Applies the same transformations as MutSig2CV:
%     * Sets null value to 200
%     * Converts null values (200) to NaN
%     * Converts chromosome names to numbers
%   - Conservation scores typically range from 0-100
%   - Higher scores indicate higher conservation
%
% Author: Generated for MutSig2CV conservation data extraction
% Date: 2024

% Input validation
if nargin ~= 5
    error('get_conservation_data requires exactly 5 input arguments');
end

if ~exist(conservation_fwb_file, 'file')
    error('Conservation FWB file not found: %s', conservation_fwb_file);
end

if ~exist(FixedWidthBinary_jar_file, 'file')
    error('FixedWidthBinary.jar file not found: %s', FixedWidthBinary_jar_file);
end

if start_pos > end_pos
    error('Start position (%d) must be <= end position (%d)', start_pos, end_pos);
end

if start_pos < 1
    error('Start position must be >= 1');
end

% Add jar to java classpath
javaclasspath(FixedWidthBinary_jar_file);

% Convert chromosome to numeric format
chr_num = convert_chr(chr);
if isnan(chr_num)
    error('Invalid chromosome identifier: %s', chr);
end

% Open FWB file
try
    fwb = org.broadinstitute.cga.tools.seq.FixedWidthBinary(conservation_fwb_file);
    fwb.setNullVal(200);
catch ME
    error('Failed to open FWB file: %s', ME.message);
end

% Extract data
try
    raw_data = double(fwb.get(chr_num, start_pos, end_pos));
catch ME
    fwb.close();
    error('Failed to extract data from FWB file: %s', ME.message);
end

% Close FWB file
fwb.close();

% Apply transformations (same as MutSig2CV)
% Convert null values (200) to NaN
raw_data(raw_data == 200) = NaN;

% Create output structure
conservation_data = struct();
conservation_data.data = raw_data;
conservation_data.positions = start_pos:end_pos;
conservation_data.chr = chr_num;
conservation_data.start = start_pos;
conservation_data.end = end_pos;
conservation_data.length = length(raw_data);

% Calculate statistics
valid_data = raw_data(~isnan(raw_data));
conservation_data.num_valid = length(valid_data);

if conservation_data.num_valid > 0
    conservation_data.mean_conservation = mean(valid_data);
    conservation_data.min_conservation = min(valid_data);
    conservation_data.max_conservation = max(valid_data);
else
    conservation_data.mean_conservation = NaN;
    conservation_data.min_conservation = NaN;
    conservation_data.max_conservation = NaN;
end

% Add summary information
conservation_data.fraction_valid = conservation_data.num_valid / conservation_data.length;
conservation_data.fraction_missing = 1 - conservation_data.fraction_valid;

end 