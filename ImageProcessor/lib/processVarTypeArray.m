function [varCellArray,varStrCellArray] = processVarTypeArray(varCellArrayStr)
% Convert variable string to the correct type

    function [var,varStr] = processVarType(inVar)
    % if the first and last characters are apostrophes ('), it is a string
    var = inVar;
    varStr = inVar;
    
    if isnumeric(inVar)
        varStr = num2str(inVar);
        return;
    end
    
    % for empty string or empty array
    if isempty(inVar)
        var = [];
        varStr = '';
        return;
    end

    if strcmp(inVar(1),'''') && strcmp(inVar(end),'''')
        var = inVar(2:end-1);
        return;
    end
    
    if strcmp(inVar(1),'{') && strcmp(inVar(end),'}')
        cellArrayStr = strsplit(inVar(2:end-1),',');
        [var,tmpStr] = processVarTypeArray(cellArrayStr);
        return;
    end

    if strcmp(inVar(1),'[') && strcmp(inVar(end),']')
        arrayStr = strsplit(inVar(2:end-1),',');
        [arrayCell,tmpStr] = processVarTypeArray(arrayStr);
        % if conversion to array is valid, return the array instead
        try
            var = cell2mat(arrayCell);
        end
        return;
    end
    
    num = str2num(inVar);
    if ~isempty(num)
        var = num;
    end
    end %function

[varCellArray,varStrCellArray] = cellfun(@processVarType,varCellArrayStr,'UniformOutput',false);

end
