% --- Executes on button press in init_pushbutton.
function handles = initDictyTracker(handles)
%%%%% Initialize variables

handles = init_guivars(handles,'DictyTracker_vars.xlsx');
handles.subVars = handles.outVars;
handles = getVars(handles,handles.subVars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables specific to this sub figure
handles.savedFields = handles.subVars.savedFields;
set(handles.plotFields_popupmenu,'String',handles.savedFields);
handles.savedFields'

% Color map for plotting multiple curves
handles.colors = lines(10);
