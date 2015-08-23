% --- Executes on button press in init_pushbutton.
function handles = init_guivars(handles,varsFileName)
%%%%% Initialize variables
handles.outVars = struct;

% Load scalar variables
vars = readVarScalars(varsFileName);
for iKey=1:length(vars.keys)
    key = vars.keys{iKey};
    handles.outVars.(key) = vars.(key);
%     updatestatus(handles.main_figure,['Set variable ' key ' specified in ' handles.varsFileName ' to: ' vars.varStr.(key)]);
end

% Load vector variables
vars = readVarVectors(handles.outVars.varVectorsFileName,1);
for iKey=1:length(vars.keys)
    key = vars.keys{iKey};
    handles.outVars.(key) = vars.(key);
%     updatestatus(handles.main_figure,['Set variable ' key ' specified in ' handles.varVectorsFileName ' to:']);
%     updatestatus(handles.main_figure,vars.varStr.(key));
end
