% Author: Duong "Zoom" Nguyen
% Email: duongnh "at" gmail "dot" com
% version 1.4

function varargout = DictyTracker(varargin)
% DICTYTRACKER MATLAB code for DictyTracker.fig
%      DICTYTRACKER, by itself, creates a new DICTYTRACKER or raises the existing
%      singleton*.
%
%      H = DICTYTRACKER returns the handle to a new DICTYTRACKER or the handle to
%      the existing singleton*.
%
%      DICTYTRACKER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICTYTRACKER.M with the given input arguments.
%
%      DICTYTRACKER('Property','Value',...) creates a new DICTYTRACKER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DictyTracker_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DictyTracker_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DictyTracker

% Last Modified by GUIDE v2.5 22-Jan-2015 13:49:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DictyTracker_OpeningFcn, ...
                   'gui_OutputFcn',  @DictyTracker_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before DictyTracker is made visible.
function DictyTracker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DictyTracker (see VARARGIN)

% Choose default command line output for ImageProcessor
handles.output = hObject;
% handles = initDictyTracker(handles);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = DictyTracker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];


% --- Executes when user attempts to close dicty_figure.
function dicty_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to dicty_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure


% uiresume(hObject);
delete(hObject);

        
% --- Executes on button press in findSpots_pushbutton.
function handles = findSpots_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to findSpots_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded to the left axes.');
    return;
end

if ~isfield(handles,'yfpFilePattern') || ~isfield(handles,'cfpFilePattern')
    updatestatus(handles,'No FRET images have been loaded.');
    return;
end

inVars = struct;
inVars.imdir = handles.filedir;
inVars.y_name = handles.yfpFilePattern;
inVars.c_name = handles.cfpFilePattern;
inVars.startImg = handles.iCurrentFile;

handles.startImg = handles.iCurrentFile;

% outVars = getFRET(inVars);
outVars = getFRETZ(inVars);

handles.tracks = outVars.tracks;
handles.masks = outVars.masks;
% handles.tottracks = outVars.tottracks;

handles.numSpots = length(handles.tracks);
if (handles.numSpots ==0) || (handles.numSpots ==1 && isempty(handles.tracks(1).x))
    updatestatus(handles,'No spot has been detected.');
    return;
end
handles.allSpotPos = [handles.tracks(:).x;handles.tracks(:).y]';

imshow(handles.imgLeft,'Parent',handles.left_axes);
handles = none_pushbutton_Callback(hObject, eventdata, handles);

if handles.isShowRightImage
    handles.imgRight = handles.imgLeft;
    allMasks = zeros(size(handles.imgRight));
    for i=1:length(handles.masks)
        allMasks = allMasks|handles.masks{i};
    end
    handles.imgRight(allMasks) = handles.maxPixelValue;

    imshow(handles.imgRight,'Parent',handles.right_axes);
    for i=1:handles.numSpots
        text(handles.allSpotPos(i,1),handles.allSpotPos(i,2),['\color{green}' num2str(i)],'Parent',handles.right_axes);
    end
    
end

guidata(hObject, handles);


% --- Executes on button press in saveTrack_pushbutton.
function saveTrack_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveTrack_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'tracks')
    updatestatus(handles,'No track data found.');
    return;
end

indices = find(handles.selectedIndexArray);

if length(indices)==0
    updatestatus(handles,'No cell has been selected for saving.');
    return;
end

updatestatus(handles,'Saving tracks data. Existing tracks will be overwritten.');

allFRET = [];
allNames = {};
for index=indices
    savedBaseName = ['track' num2str(index)];
    savedFilePath = [handles.filedir savedBaseName '.txt'];
    allNames{end+1} = savedBaseName;
    
    % Save only the first selected track
    trackToSave = handles.tracks(index);
    
    allFRET = [allFRET trackToSave.ratio(:)];
    
    trackToSaveData = [];
    for iField=1:length(handles.savedFields)
        field = handles.savedFields{iField};
        fielddata = trackToSave.(field);
        trackToSaveData = [trackToSaveData fielddata(:)];
    end

    saveStringArray(savedFilePath,handles.savedFields);
    dlmwrite(savedFilePath, trackToSaveData,'delimiter', '\t', '-append');
    
    posStr = ['(' num2str(round(handles.allSpotPos(index,1))) ',' num2str(round(handles.allSpotPos(index,2))) ')'];
    updatestatus(handles,['Track number ' num2str(index) ' with starting position ' posStr ' has been saved.']);
end

allPath = fullfile(handles.filedir,'AllFRET.txt');
saveStringArray(allPath,allNames);
dlmwrite(allPath, allFRET,'delimiter', '\t','-append');

guidata(hObject,handles)


function pix2um_edit_Callback(hObject, eventdata, handles)
% hObject    handle to pix2um_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of pix2um_edit as text
%        str2double(get(hObject,'String')) returns contents of pix2um_edit as a double
handles.pix2um = str2double(get(hObject,'String'));
updatestatus(handles,['Calibration is set to: 1 pixel = ' num2str(handles.pix2um) 'um']);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pix2um_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pix2um_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function cfpSuffix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to cfpSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of cfpSuffix_edit as text
%        str2double(get(hObject,'String')) returns contents of cfpSuffix_edit as a double
handles.cfpSuffix = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function cfpSuffix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cfpSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function handles = TrackSelectedSpots_Zoom(handles,inVars)
% Track cells by choosing ones that are closest to the cells in the
% previous frame

% Status message to return
status = '';

outTracks = inVars.tracks;

% Setting up the waitbar.
w = waitbar(0, 'Please wait...');
message = 'Starting ...';
waitbar(0, w, message);

% Positions of the cells that have been identified and selected
currentSpotPos = [inVars.tracks(:).x;inVars.tracks(:).y]';

startFrame = handles.iCurrentFile;
while handles.iCurrentFile < handles.numFiles
    % if waitbar is closed, cancel the loop
    if(~ishandle(w))
        status = 'Waitbar is closed. Tracking is halted.';
        break;
    end
    
    % Displaying the waitbar.
    message = ['Analyzing Frame ' num2str(handles.iCurrentFile) '. Last frame to analyze: ' num2str(handles.numFiles) '.'];
    waitbar((handles.iCurrentFile-startFrame)/(handles.numFiles-startFrame), w, message);

    handles.iCurrentFile = handles.iCurrentFile+1;
    updatestatus(handles,['Processing frame ' num2str(handles.iCurrentFile)]);
    
    set(handles.current_edit,'String',num2str(handles.iCurrentFile));
    handles = handles.current_edit_Callback(handles.current_edit, 0, handles);
    
    handles = findSpots_pushbutton_Callback(handles.findSpots_pushbutton, 0, handles);
    hold(handles.left_axes,'on');
    for iSpot=1:size(currentSpotPos,1)
        iNear = nearest(handles.allSpotPos,currentSpotPos(iSpot,:));
        % add the values in handles.tracks to outTracks
        allFields = fields(outTracks);
        for i=1:length(allFields)
            f = allFields{i};
            tmp = outTracks(iSpot).(f);
            tmp(end+1)=handles.tracks(iNear).(f);
            outTracks(iSpot).(f) = tmp;
        end
        currentSpotPos(iSpot,:) = handles.allSpotPos(iNear,:);
        plot(handles.left_axes,currentSpotPos(iSpot,1),currentSpotPos(iSpot,2),'g+');
    end
    hold(handles.left_axes,'off');
    
end
handles.tracks = outTracks;

if ishandle(w)
    close(w); %waitbar
end

status = [status 'Tracking has finished.'];
setappdata(handles.main_figure,'status',status);


function handles = gridFRET(handles,inVars)
% Divide frame into divisions and track FRET signal on each over time

% Status message to return
status = '';

outTracks = repmat(struct,handles.numDivisions^2,1);

% Setting up the waitbar.
w = waitbar(0, 'Please wait...');
message = 'Starting ...';
waitbar(0, w, message);

startFrame = handles.iCurrentFile;
while handles.iCurrentFile <= handles.numFiles
    % if waitbar is closed, cancel the loop
    if(~ishandle(w))
        status = 'Waitbar is closed. Tracking is halted.';
        break;
    end
    
    % Displaying the waitbar.
    message = ['Analyzing Frame ' num2str(handles.iCurrentFile) '. Last frame to analyze: ' num2str(handles.numFiles) '.'];
    waitbar((handles.iCurrentFile-startFrame)/(handles.numFiles-startFrame), w, message);

    updatestatus(handles,['Processing frame ' num2str(handles.iCurrentFile)]);
    
    set(handles.current_edit,'String',num2str(handles.iCurrentFile));
    handles = handles.current_edit_Callback(handles.current_edit, 0, handles);
      
    inVars.startImg = handles.iCurrentFile;
    output_single = gridFRET_single(inVars);
    
    handles.imgRight = handles.imgLeft;
    handles.imgRight(output_single.bw) = handles.maxPixelValue;
    imshow(handles.imgRight,'Parent',handles.right_axes);

    for iSpot=1:handles.numDivisions^2
        % add the values in handles.tracks to outTracks
        allFields = fields(output_single.tracks);
        for i=1:length(allFields)
            f = allFields{i};
            if ~isfield(outTracks(iSpot),f)
                outTracks(iSpot).(f) = output_single.tracks(iSpot).(f);
            else
                outTracks(iSpot).(f) = [outTracks(iSpot).(f) output_single.tracks(iSpot).(f)];
            end
        end
    end
    
    if inVars.isSaveFRETImages
        outfilePath = fullfile(handles.filedir,['img_' num2str(handles.iCurrentFile,'%09d') '_FRET.tif']);
%         tmpimg = autocontrast(output_single.Iratio,16);
        imwrite(output_single.Iratio,outfilePath);
    end
    
    handles.iCurrentFile = handles.iCurrentFile+1;
end
handles.tracks = outTracks;

if ishandle(w)
    close(w); %waitbar
end
 
status = [status 'Tracking has finished.'];
setappdata(handles.main_figure,'status',status);


% --- Executes on button press in track_pushbutton.
function handles = track_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to track_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded.');
    return;
end

if ~isfield(handles,'numSpots')
    updatestatus(handles,'No cell has been detected.');
    return;
end

% indices of cells to remove
removedIndices = find(1-handles.selectedIndexArray);
trackedIndices = find(handles.selectedIndexArray);

if length(removedIndices)>=handles.numSpots
    updatestatus(handles,'No cell has been selected for tracking.');
    return;
end

inVars = struct;
% inVars.tottracks = handles.tottracks;
inVars.tracks = handles.tracks;
inVars.masks = handles.masks;
inVars.tracks(removedIndices) = [];
inVars.masks(removedIndices) = [];

inVars.dt = 1;
inVars.tau = 1;
inVars.tau2 = 1;
inVars.imdir = handles.filedir;
inVars.y_name = handles.yfpFilePattern;
inVars.c_name = handles.cfpFilePattern;
inVars.startImg = handles.startImg;
inVars.numImg = handles.numFiles-handles.startImg+1;

% Darvin's tracking
% outVars = TrackSelectedSpots(inVars);

handles = TrackSelectedSpots_Zoom(handles,inVars);
updatestatus(handles);

handles.numSpots = length(handles.tracks);

% Initial position of each track
handles.allSpotPos = [];
for index=1:length(handles.tracks)
    track = handles.tracks(index);
    handles.allSpotPos = [handles.allSpotPos;[track.x(1) track.y(1)]];
    posStr = ['(' num2str(round(track.x(1))) ',' num2str(round(track.y(1))) ')'];
    updatestatus(handles,['Spot number ' num2str(trackedIndices(index)) ' with starting position ' posStr ' has been tracked.']);
end

% Mark all tracked cells
handles = all_pushbutton_Callback(hObject, eventdata, handles);
drawTracks_pushbutton_Callback(hObject, eventdata, handles);

if handles.isAutoSaveTrack
    saveTrack_pushbutton_Callback(hObject, eventdata, handles);
end

guidata(hObject, handles);


% --- Executes on button press in all_pushbutton.
function handles = all_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to all_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'numSpots')
    updatestatus(handles,'No cell has been detected.');
    return;
end

% array to keep track of which cells are selected with value of 1. Zero
% means not selected
handles.selectedIndexArray = ones(1,handles.numSpots);

hold(handles.left_axes,'on');
plot(handles.left_axes,handles.allSpotPos(:,1),handles.allSpotPos(:,2),'g+');
hold(handles.left_axes,'off');

guidata(hObject, handles);


% --- Executes on button press in none_pushbutton.
function handles = none_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to none_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'numSpots')
    updatestatus(handles,'No cell has been detected.');
    return;
end

% array to keep track of which cells are selected with value of 1. Zero
% means not selected
handles.selectedIndexArray = zeros(1,handles.numSpots);

hold(handles.left_axes,'on');
plot(handles.left_axes,handles.allSpotPos(:,1),handles.allSpotPos(:,2),'r+');
hold(handles.left_axes,'off');

guidata(hObject, handles);


% --- Executes on button press in drawTracks_pushbutton.
function drawTracks_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to drawTracks_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
indices = find(handles.selectedIndexArray);
if length(indices)==0
    updatestatus(handles,'No cell has been selected.');
    return;
end

% if there is only one track, color it with jet colormap
% if length(indices)==1
%     track = [handles.tracks(indices).x(:) handles.tracks(indices).y(:)];
%     numPoints = length(track(:,1));
%     colors = jet(handles.numPoints);
% 
%     for i = 1:numPoints-1
%         color = colors(mod(i-1,handles.numPoints)+1,:);
%         line(handles.left_axes,track(i:i+1,1),track(i:i+1,2),'Color',color,'LineWidth',3);
%     end
% else

%    axes(handles.left_axes);

% invert image to that background is white to save ink!
handles.imgLeft = handles.maxPixelValue-handles.imgLeft;

imshow(handles.imgLeft,'Parent',handles.left_axes);
hold(handles.left_axes,'on');

% plot(handles.left_axes,handles.allSpotPos(:,1),handles.allSpotPos(:,2),'r+');

numColors = length(handles.colors);
for iTrack = indices
    color = handles.colors(mod(iTrack-1,numColors)+1,:);
    plot(handles.left_axes,handles.tracks(iTrack).x,handles.tracks(iTrack).y,'-','Color',color,'LineWidth',2);
%     plot(handles.left_axes,handles.allSpotPos(iTrack,1),handles.allSpotPos(iTrack,2),'g+');
    text(handles.allSpotPos(iTrack,1),handles.allSpotPos(iTrack,2),['\color{cyan}' num2str(iTrack)],'Parent',handles.left_axes);
end

hold(handles.left_axes,'off');



% --- Executes on selection change in plotFields_popupmenu.
function plotFields_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotFields_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotFields_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotFields_popupmenu
handles.iPlotField = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function plotFields_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotFields_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in plot_pushbutton.
function plot_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plot_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'numSpots')
    updatestatus(handles,'No cell has been detected.');
    return;
end

indices = find(handles.selectedIndexArray);
if length(indices)==0
    updatestatus(handles,'No cell has been selected.');
    return;
end
plotField = handles.savedFields{handles.iPlotField};

cla(handles.right_axes,'reset');
hold(handles.right_axes,'on');

xlabel(handles.right_axes,'time (min)');
ylabel(handles.right_axes,plotField);

numColors = length(handles.colors);
for iTrack = indices
    color = handles.colors(mod(iTrack-1,numColors)+1,:);
    time = handles.tracks(iTrack).frame./handles.framePerMin;
    plot(handles.right_axes,time,handles.tracks(iTrack).(plotField),'-','Color',color,'LineWidth',2);
end

hold(handles.right_axes,'off');


function framePerMin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to framePerMin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of framePerMin_edit as text
%        str2double(get(hObject,'String')) returns contents of framePerMin_edit as a double
handles.framePerMin = str2double(get(hObject,'String'));
updatestatus(handles,['Frame rate is set to: ' num2str(handles.framePerMin) ' frame(s) = 1 min']);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function framePerMin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to framePerMin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in grid_pushbutton.
function handles = grid_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to grid_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded.');
    return;
end

handles.startImg = handles.iCurrentFile;

inVars.imdir = handles.filedir;
inVars.imdir
inVars.y_name = handles.yfpFilePattern;
inVars.c_name = handles.cfpFilePattern;
inVars.threshold_name = handles.yfpFilePattern;
% inVars.startImg = handles.startImg;
inVars.numImg = handles.numFiles-handles.startImg+1;
inVars.numRowDiv = handles.numDivisions;
inVars.numColDiv = handles.numDivisions;

if handles.numDivisions==1
    inVars.isSaveFRETImages = true;
else 
    inVars.isSaveFRETImages = false;
end

% Do not threshold before taking the ratio c/y
inVars.isThresholdFRET = handles.isThresholdFRET;
inVars.isThresholdAdaptive = handles.isThresholdAdaptive;
inVars.threshold = handles.threshold;

handles = gridFRET(handles,inVars);
updatestatus(handles);

% The following variables are set for compatibility with regular tracking
handles.numSpots = length(handles.tracks);
handles.selectedIndexArray = ones(1,handles.numSpots);

handles.allSpotPos = [];
for index=1:length(handles.tracks)
    track = handles.tracks(index);
    handles.allSpotPos = [handles.allSpotPos;[track.x(1) track.y(1)]];
end

if handles.isAutoSaveTrack
    saveTrack_pushbutton_Callback(hObject, eventdata, handles);
end

guidata(hObject, handles);

function numDivisions_edit_Callback(hObject, eventdata, handles)
% hObject    handle to numDivisions_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numDivisions_edit as text
%        str2double(get(hObject,'String')) returns contents of numDivisions_edit as a double
handles.numDivisions = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function numDivisions_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numDivisions_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function tolerance_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tolerance_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tolerance_edit as text
%        str2double(get(hObject,'String')) returns contents of tolerance_edit as a double
handles.tolerance = str2double(get(hObject,'String'));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tolerance_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tolerance_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in magicWand_pushbutton.
function magicWand_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to magicWand_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
bin_mask = magicwand_gray(handles.imgLeft, handles.selectedPosArray(:,2), handles.selectedPosArray(:,1), handles.tolerance);
handles.imgRight = handles.imgLeft;
handles.imgRight(bin_mask) = handles.maxPixelValue;
imshow(handles.imgRight,'Parent',handles.right_axes);


% --- Executes on button press in loadFRET_pushbutton.
function handles = loadFRET_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadFRET_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in load_pushbutton.
handles.yfpFilePattern = [handles.prefix '*' handles.yfpSuffix '.' handles.imgExtension];
handles.cfpFilePattern = [handles.prefix '*' handles.cfpSuffix '.' handles.imgExtension];
% handles.fretFilePattern = [handles.prefix '*' handles.fretSuffix '.' handles.imgExtension];

outVars = handles.loadFiles(handles,handles.filedir,handles.cfpFilePattern);
handles.cfpFiles = outVars.files;
handles.cfpFileNames = outVars.fileNames;
handles.cfpNumFiles = outVars.numFiles;

% if ~isempty(handles.fretSuffix)
%     outVars = handles.loadFiles(handles,handles.filedir,handles.fretFilePattern);
%     handles.fretFiles = outVars.files;
%     handles.fretFileNames = outVars.fileNames;
%     handles.fretNumFiles = outVars.numFiles;
% end

handles = handles.load_pushbutton_Callback(hObject, eventdata, handles);
if handles.numFiles ~= handles.cfpNumFiles
    updatestatus(handles,'Unequal number of YFP and CFP files.');
    return;
end

guidata(hObject, handles);


% --- Executes on button press in isAutoSaveTrack_checkbox.
function isAutoSaveTrack_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isAutoSaveTrack_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isAutoSaveTrack_checkbox
handles.isAutoSaveTrack = get(hObject,'Value');
guidata(hObject, handles);



function yfpSuffix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to yfpSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yfpSuffix_edit as text
%        str2double(get(hObject,'String')) returns contents of yfpSuffix_edit as a double
handles.yfpSuffix = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yfpSuffix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yfpSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in mainToYFP_pushbutton.
function mainToYFP_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to mainToYFP_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.yfpSuffix = handles.suffix;
set(handles.yfpSuffix_edit,'String',handles.yfpSuffix);
guidata(hObject, handles);


% --- Executes on button press in yfpToCFP_pushbutton.
function yfpToCFP_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to yfpToCFP_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.cfpSuffix = handles.yfpSuffix;
set(handles.cfpSuffix_edit,'String',handles.cfpSuffix);
% handles.fretSuffix = handles.yfpSuffix;
% set(handles.fretSuffix_edit,'String',handles.fretSuffix);
guidata(hObject, handles);


% --- Executes on button press in isThresholdFRET_checkbox.
function isThresholdFRET_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isThresholdFRET_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isThresholdFRET_checkbox
handles.isThresholdFRET = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in isThresholdAdaptive_checkbox.
function isThresholdAdaptive_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isThresholdAdaptive_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isThresholdAdaptive_checkbox
handles.isThresholdAdaptive = get(hObject,'Value');
guidata(hObject, handles);


function threshold_edit_Callback(hObject, eventdata, handles)
% hObject    handle to threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold_edit as text
%        str2double(get(hObject,'String')) returns contents of threshold_edit as a double
handles.threshold = str2double(get(hObject,'String'))/handles.maxPixelValue;
tmpimg = im2bw(handles.imgLeft,handles.threshold);
handles.imgRight = bwareaopen(tmpimg, 4, 4);

imshow(handles.imgRight,'Parent',handles.right_axes);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function threshold_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function fretSuffix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fretSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fretSuffix_edit as text
%        str2double(get(hObject,'String')) returns contents of fretSuffix_edit as a double


% --- Executes during object creation, after setting all properties.
function fretSuffix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fretSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in processQueue_pushbutton.
function processQueue_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to processQueue_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if length(handles.dirQueue)<1
    updatestatus(handles,'Queue is empty.');
    return;
end

filedir_orig = handles.filedir;
for iDir = 1:length(handles.dirQueue)
    handles.filedir = handles.dirQueue{iDir};
    set(handles.filedir_edit,'String',handles.filedir);
    handles = loadFRET_pushbutton_Callback(hObject, eventdata, handles);
    pause(1);
    handles = grid_pushbutton_Callback(hObject, eventdata, handles);
end

handles.filedir = filedir_orig;
set(handles.filedir_edit,'String',handles.filedir);
handles.dirQueue = {};

guidata(hObject, handles);
