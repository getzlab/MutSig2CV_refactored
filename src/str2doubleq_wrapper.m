function out = str2doubleq_wrapper(in)
% STR2DOUBLEQ_WRAPPER - Pure MATLAB implementation of str2doubleq
% This replaces the MEX-based str2doubleq function for compatibility

if iscell(in)
    out = cell(size(in));
    for i = 1:numel(in)
        if ischar(in{i})
            out{i} = str2double(in{i});
        else
            out{i} = NaN;
        end
    end
    out = cell2mat(out);
elseif ischar(in)
    out = str2double(in);
else
    out = NaN;
end

% Handle complex numbers if present
if ~isreal(out)
    out = real(out);
end
