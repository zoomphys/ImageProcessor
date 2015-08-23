function vars = readVarVectors(filepathVarVectors,isColHeaders)
% Read vector variables from filepathVarVectors. Each column contains the
% values for each variable with the name specified in the first row if
% isColHeaders is true
[nData, tData, allData] = xlsread(filepathVarVectors);
[numRows,numCols] = size(allData);

if isColHeaders
    varNames = allData(1,:);
    startRow = 2;
else
    varNames = arrayfun(@(x) ['var' num2str(x)],1:numCols,'UniformOutput',false);
    startRow = 1;
end

if startRow>numRows
    vars = [];
    return;
end

vars = struct;
vars.varStr = struct;
vars.keys = varNames;

% Populate variable with vectors. All empty elements are discarded
for iVar=1:numCols
    varName = varNames{iVar};
    varVector = allData(startRow:end,iVar);
    fh = @(x) all(isnan(x(:)));
    varVector(cellfun(fh,varVector)) = [];
    [varVector,varStrVector] = processVarTypeArray(varVector');
    vars.(varName) = varVector;
    vars.varStr.(varName) = strjoin(varStrVector,',');
end
end
