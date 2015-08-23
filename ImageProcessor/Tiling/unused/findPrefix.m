function [prefix,suffix] = findPrefix(filename,pattern)
%Find the position of the first numeric char in a consecutive string of numeric chars at the end. Return the prefix before the found position and the string containing the number.

    %position of the first non numeric char in the reverse string
    if isempty(pattern)
        pattern = '[^0-9]';
    end
    
    posFromEnd = regexp(fliplr(filename),pattern,'once');
    if isempty(posFromEnd)
        posFromStart = 1;
    else
        posFromStart = length(filename)-posFromEnd+2;
    end

    prefix = filename(1:posFromStart-1);
    suffix = filename(posFromStart:end);
end

