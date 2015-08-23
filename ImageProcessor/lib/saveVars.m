function handles = saveVars(handles,varStruct)
% save gui variables
guiVars = struct;
for i=1:length(varStruct.savedVars)
    var = varStruct.savedVars{i};
    if isfield(handles,var)
        guiVars.(var) = handles.(var);
    end
end
save(varStruct.guivarsFileName,'guiVars');
