function vars = readVarScalars(filepathVarVectors)
% Read vector variables from filepathVarVectors. Each column contains the
% values for each variable with the name specified in the first row if
% isColHeaders is true
[nData, tData, allData] = xlsread(filepathVarVectors);
[numRows,numCols] = size(allData);

vars = struct;
vars.varStr = struct; % store string representation of variable
vars.keys = allData(:,1)';

varVector = allData(:,2)';
[varVector,varStrVector] = processVarTypeArray(varVector);

% Populate variable with values. All empty elements are discarded
for iVar=1:numRows
    varName = vars.keys{iVar};
    vars.(varName) = varVector{iVar};
    vars.varStr.(varName) = varStrVector{iVar};
end
