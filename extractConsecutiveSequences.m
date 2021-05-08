function [out cat_out] = extractConsecutiveSequences(in, gap)
% keep consecutive values with tolerance = gap
%
% Input:
% ------
% in            row vector
% gap           maximum allowed gap in days between two values
%
% Output:
% -------
% out           cell array with one series per column
% cat_out       concatenated series

    % find where adjacent elements are consecutive
    id = [false, diff(in)<=gap+1, false];
    % find idx of beginning of consecutive numbers
    start = strfind(id,[0 1]);
    % find idx of end
    stop = strfind(id,[1 0]);
    
    out = cell(1,length(start));
    cat_out = [];
    % write each sequence separately
    for i = 1:length(start) 
        out{i} = in(start(i):stop(i));
        cat_out = cat(2,cat_out,in(start(i):stop(i)));
    end
    cat_out = unique(cat_out);
end