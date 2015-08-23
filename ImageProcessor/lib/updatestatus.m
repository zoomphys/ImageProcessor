function updatestatus(varargin)
% update the status bar in the GUI

if nargin==0 || nargin>2
    disp('The number of arguments need to be 1 or 2.');
end

handles = varargin{1};

if nargin==1
    msg = getappdata(handles.main_figure,'status');
    setappdata(handles.main_figure,'status','');
else
    msg = varargin{2};
end

if isempty(msg)
    msg = 'No status to display.';
end

log = getappdata(handles.main_figure,'log');
log = {log{:} [getTime() '   ' msg]};
% log = {log{:} msg};

disp(msg);
%msg
set(handles.status_text,'String',msg);
setappdata(handles.main_figure,'log',log);
