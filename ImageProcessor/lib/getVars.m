function handles = getVars(handles,varStruct)
% get gui variables from a file specified by varStruct and set the
% corresponding gui elements

%Only to bootstrap guiVars
%guiVars = struct;

% Initialize variables from guiVars file
if exist(varStruct.guivarsFileName ,'file')
    load(varStruct.guivarsFileName ,'guiVars');
end
for i=1:length(varStruct.savedVars)
    var = varStruct.savedVars{i};
    if isfield(guiVars,var)
        handles.(var) = guiVars.(var);
        guiName = varStruct.guiName{i};
        if ~isempty(guiName)
            set(handles.(guiName),varStruct.guiType{i},handles.(var));
        end
    end
end
