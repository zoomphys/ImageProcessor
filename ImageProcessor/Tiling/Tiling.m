% Author: Duong "Zoom" Nguyen
% Email: duongnh "at" gmail "dot" com
% version 1.4

function varargout = Tiling(varargin)
% TILING MATLAB code for Tiling.fig
%      TILING, by itself, creates a new TILING or raises the existing
%      singleton*.
%
%      H = TILING returns the handle to a new TILING or the handle to
%      the existing singleton*.
%
%      TILING('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TILING.M with the given input arguments.
%
%      TILING('Property','Value',...) creates a new TILING or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Tiling_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Tiling_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Tiling

% Last Modified by GUIDE v2.5 29-Jan-2015 10:28:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Tiling_OpeningFcn, ...
                   'gui_OutputFcn',  @Tiling_OutputFcn, ...
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


% --- Executes just before Tiling is made visible.
function Tiling_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Tiling (see VARARGIN)

% Choose default command line output for ImageProcessor
handles.output = hObject;
% handles = initDictyTracker(handles);

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = Tiling_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% varargout{1} = handles.output;
varargout{1} = [];


% --- Executes when user attempts to close subGUI_figure.
function subGUI_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to subGUI_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure


% uiresume(hObject);
delete(hObject);
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Executes on button press in replMainSuffix_pushbutton.
function replMainSuffix_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to replMainSuffix_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% handles.suffix = [handles.profiles{handles.iProfile} '_' handles.tilingSuffix];
handles.suffix = handles.profiles{handles.iProfile};
set(handles.suffix_edit,'String',handles.suffix);
guidata(hObject, handles);



function tilingSuffix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to tilingSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tilingSuffix_edit as text
%        str2double(get(hObject,'String')) returns contents of tilingSuffix_edit as a double
handles.tilingSuffix = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function tilingSuffix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tilingSuffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in stitch_pushbutton.
function stitch_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to stitch_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded.');
    return;
end
status = '';

% filedir = 'E:\Archive\AllRawData\2013-07-12_LaserUncaging_d-cA1uMCellsChemotax\set2\image_1';
trackfile = 'MMControl_track.txt';
% pix2um = 0.162;
pix2um = 0.162;

startFrame = handles.iCurrentFile;
track = importdata(fullfile(handles.filedir,trackfile));
xShift = track.data(:,1)/pix2um;
yShift = track.data(:,2)/pix2um;
xCell = track.data(:,3)/pix2um;
yCell = track.data(:,4)/pix2um;

nTrack = length(xShift);

x = zeros(nTrack,1);
y = zeros(nTrack,1);
for i=2:nTrack
    x(i) = x(i-1)+xShift(i-1);
    y(i) = y(i-1)+yShift(i-1);
end
xCell = xCell+x;
yCell = yCell+y;
trajTotal = [xCell(:) yCell(:)]

xmin = floor(min(x));
xmax = ceil(max(x));
ymin = floor(min(y));
ymax = ceil(max(y));

profileList = {handles.suffix};
% profileList = {'CFP2-YFP_000','CFP2-CFP_000','GFP_000','Transmission_000'};

% Setting up the waitbar.
w = waitbar(0, 'Please wait...');
message = 'Starting ...';
waitbar(0, w, message);

% while handles.iCurrentFile <= startFrame+1
while handles.iCurrentFile <= handles.numFiles
    % if waitbar is closed, cancel the loop
    if(~ishandle(w))
        status = 'Waitbar is closed. Stitching is halted.';
        break;
    end
    
    % Displaying the waitbar.
    message = ['Analyzing Frame ' num2str(handles.iCurrentFile) '. Last frame to analyze: ' num2str(handles.numFiles) '.'];
    waitbar((handles.iCurrentFile-startFrame)/(handles.numFiles-startFrame), w, message);

    set(handles.current_edit,'String',num2str(handles.iCurrentFile));
    handles = handles.current_edit_Callback(handles.current_edit, 0, handles);
    updatestatus(handles,['Processing frame ' num2str(handles.iCurrentFile)]);
    
    % translation in units of pixels
    
    dx = floor(x(handles.iCurrentFile));
    dy = floor(y(handles.iCurrentFile));
    traj = trajTotal(startFrame:handles.iCurrentFile,:)
    
    for profile = profileList
        profileStr = profile{1};
        filename = ['img_' num2str(handles.iCurrentFile,'%09d') '_' profileStr '.tif'];
        
        im = imread(fullfile(handles.filedir,filename));
        [nRows, nCols] = size(im);
        bw = im2bw(im, graythresh(im));
        bkgd = floor(mean(double(im(~bw))));
        
        tmpimg = uint16(ones(nRows-ymin+ymax,nCols-xmin+xmax))*bkgd;
        tmpimg(1+dy-ymin:nRows+dy-ymin,1+dx-xmin:nCols+dx-xmin) = im;
        
        handles.imgRight = tmpimg;
        if handles.isAutoContrast
            handles.imgRight = autocontrast(handles.imgRight,16);
        end
        imshow(handles.imgRight,'Parent',handles.right_axes);

        % draw track
%         numPoints = size(traj,1);
%         colors = jet(numPoints);
% 
%         if numPoints>1
%             for i = 1:numPoints-1
%                 color = colors(i,:);
%                 line(traj(i:i+1,1),traj(i:i+1,2),'Color',color,'LineWidth',3,'Parent',handles.right_axes);
%             end
%         end
        if size(traj,1)>1;
            hold(handles.right_axes,'on');
            plot(handles.right_axes,traj(:,1),traj(:,2),'y-','LineWidth',2)
            tmpimg = getframe(handles.right_axes);
            hold(handles.right_axes,'off');
        end
        outfilePath = fullfile(handles.filedir,['stitched_' num2str(handles.iCurrentFile,'%09d') '_' profileStr '.tif']);
        imwrite(tmpimg,outfilePath);
    end
    
    handles.iCurrentFile = handles.iCurrentFile+1;
end

if ishandle(w)
    close(w); %waitbar
end
 
status = [status 'Stitching has finished.'];
updatestatus(handles,status);

guidata(hObject, handles);


% --- Executes on button press in tile_pushbutton.
function tile_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tile_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded.');
    return;
end
status = '';

numSubdirs = length(handles.dirQueue);
handles.filedirTiled = fullfile(handles.maindir,'tiled');

profileList = {handles.profiles{handles.iProfile}};
if ~exist(handles.filedirTiled,'dir')
    status = mkdir(handles.filedirTiled);
    updatestatus(handles,'Output directory has been created.');
end

% Setting up the waitbar.
w = waitbar(0, 'Please wait...');
message = 'Starting ...';
waitbar(0, w, message);

startFrame = handles.iCurrentFile;
% while handles.iCurrentFile <= startFrame+1
while handles.iCurrentFile <= handles.numFiles
% while handles.iCurrentFile <= 123
    % if waitbar is closed, cancel the loop
    if(~ishandle(w))
        status = 'Waitbar is closed. Tiling is halted.';
        break;
    end
    
    % Displaying the waitbar.
    message = ['Analyzing Frame ' num2str(handles.iCurrentFile) '. Last frame to analyze: ' num2str(handles.numFiles) '.'];
    waitbar((handles.iCurrentFile-startFrame)/(handles.numFiles-startFrame), w, message);

    set(handles.current_edit,'String',num2str(handles.iCurrentFile));
    handles = handles.current_edit_Callback(handles.current_edit, 0, handles);
    updatestatus(handles,['Processing frame ' num2str(handles.iCurrentFile)]);
    
    for profile = profileList
        profileStr = profile{1};
        filename = ['img_' num2str(handles.iCurrentFile,'%09d') '_' profileStr '.tif'];
        
        imgs = {};
        for iSubdir = 1:numSubdirs
            imgs{end+1} = imread(fullfile(handles.dirQueue{iSubdir},filename));
        end
        
        imgTiled = []
        for iDir = 1:numSubdirs
%             imgTiled = [imgs{iDir},imgTiled];
            imgTiled = [imgTiled,imgs{iDir}];
        end

%         handles.imgRight = imgTiled;
%         if handles.isAutoContrast
%             handles.imgRight = autocontrast(handles.imgRight,16);
%         end
%         imshow(handles.imgRight,'Parent',handles.right_axes);

        outfilePath = fullfile(handles.filedirTiled,['img_' num2str(handles.iCurrentFile,'%09d') '_' profileStr '_' handles.tilingSuffix '.tif']);
        imwrite(imgTiled,outfilePath);
    end
    
    handles.iCurrentFile = handles.iCurrentFile+1;
end

if ishandle(w)
    close(w); %waitbar
end
 
status = [status 'Tiling has finished.'];
updatestatus(handles,status);

guidata(hObject, handles);


% --- Executes on button press in openInFiji_pushbutton.
function openInFiji_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to openInFiji_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
profileName = handles.profiles{handles.iProfile};
filePattern = [handles.prefix '*' profileName '_' handles.tilingSuffix '.' handles.imgExtension];
files = dir(fullfile(handles.filedir,filePattern));
filepath = strrep(fullfile(handles.filedir,files(1).name),'\','\\');
updatestatus(handles,['Read image sequence in Fiji starting with: ' filepath]);

repl_array_fiji = {'$filepath',filepath;
    '$numfiles',num2str(length(files));
    '$pattern',profileName};

[filedir,filename,ext] = fileparts(mfilename('fullpath'));
filepath_FijiOpenImageSequenceTemplate = fullfile(filedir,'FijiOpenImageSequenceTemplate.ijm');
filepath_FijiOpenImageSequence = fullfile(handles.filedir,'FijiOpenImageSequence.ijm');
replace_in_file(filepath_FijiOpenImageSequenceTemplate, filepath_FijiOpenImageSequence, repl_array_fiji);

system(['C:\MyPrograms\Fiji.app\ImageJ-win32.exe -macro "' filepath_FijiOpenImageSequence '"']);

repl_array_annotate = {'$umScale',handles.umScale;
    '$frmPerMin',handles.frmPerMin};

filepath_FijiAnnotateColorTemplate = fullfile(filedir,'FijiAnnotateColorTemplate.ijm');
filepath_FijiAnnotateColor = fullfile(handles.filedir,'FijiAnnotateColor.ijm');
replace_in_file(filepath_FijiAnnotateColorTemplate, filepath_FijiAnnotateColor, repl_array_annotate);


% --- Executes on button press in toMovie_pushbutton.
function toMovie_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to toMovie_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filedir_pwd = pwd;
cd(handles.filedir);
system(['pic2vid p%09d.tif ' handles.vidFPS]);
cd(filedir_pwd);

if handles.isAutoClean
    cleanOutputs_pushbutton_Callback(hObject, eventdata, handles);
end


function umScale_edit_Callback(hObject, eventdata, handles)
% hObject    handle to umScale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of umScale_edit as text
%        str2double(get(hObject,'String')) returns contents of umScale_edit as a double
handles.umScale = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function umScale_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to umScale_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function frmPerMin_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frmPerMin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frmPerMin_edit as text
%        str2double(get(hObject,'String')) returns contents of frmPerMin_edit as a double
handles.frmPerMin = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function frmPerMin_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frmPerMin_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function vidFPS_edit_Callback(hObject, eventdata, handles)
% hObject    handle to vidFPS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vidFPS_edit as text
%        str2double(get(hObject,'String')) returns contents of vidFPS_edit as a double
handles.vidFPS = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function vidFPS_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vidFPS_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cleanOutputs_pushbutton.
function cleanOutputs_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cleanOutputs_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filedir_pwd = pwd;
cd(handles.filedir);
system('del p000??????.tif');
cd(filedir_pwd);


% --- Executes on button press in tiledDirToMain_pushbutton.
function tiledDirToMain_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tiledDirToMain_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'filedirTiled') && exist(handles.filedirTiled,'dir')
    handles.filedir = handles.filedirTiled;
    set(handles.filedir_edit,'String',handles.filedir);
end
guidata(hObject, handles);


% --- Executes on selection change in profiles_popupmenu.
function profiles_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to profiles_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns profiles_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from profiles_popupmenu
handles.iProfile = get(hObject,'Value');
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function profiles_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to profiles_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadTiled_pushbutton.
function loadTiled_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadTiled_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.maindir = handles.filedir;
handles.filedir = handles.dirQueue{1};
set(handles.filedir_edit,'String',handles.filedir);

handles = handles.load_pushbutton_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


% --- Executes on button press in isAutoClean_checkbox.
function isAutoClean_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isAutoClean_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isAutoClean_checkbox
handles.isAutoClean = get(hObject,'Value');
guidata(hObject, handles);


% --- Executes on button press in projectImg_pushbutton.
function projectImg_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to projectImg_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Setting up the waitbar.
if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded.');
    return;
end
status = '';

w = waitbar(0, 'Please wait...');
message = 'Starting ...';
waitbar(0, w, message);

if handles.isCrop
    yBottom = handles.yBottomCrop;
    yTop = handles.yTopCrop;
else
    yBottom = size(handles.imgLeft,1);
    yTop = 1
end
numRows = yBottom-yTop+1;
imgProjected = [];

startFrame = handles.iCurrentFile;
while handles.iCurrentFile <= handles.numFiles
% while handles.iCurrentFile <= 9
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
    
    % only find bright pixels belonging to a cell (not background) and take
    % the average of those along the same vertical line
    img = imread(fullfile(handles.filedir,handles.fileNames{handles.iCurrentFile}));
    bw = im2bw(img, graythresh(img));
    % number of bright pixels on each vertical line, needs to be at least 1
    numPixels = max(1,sum(bw(yTop:yBottom,:)));
    line = sum(img(yTop:yBottom,:))./numPixels;
    imgProjected = [imgProjected;line];
    handles.iCurrentFile = handles.iCurrentFile+1;
end

save(fullfile(handles.filedir,'imgProjected.tif'),'imgProjected');
dlmwrite(fullfile(handles.filedir,'imgProjected.txt'),imgProjected,'\t');
hdf5write(fullfile(handles.filedir,'imgProjected.dat'), '/main', imgProjected);

imgProjected = uint16(round(imgProjected'));
imshow(imgProjected,'Parent',handles.right_axes);
imwrite(imgProjected,fullfile(handles.filedir,'imgProjected.tif'));


if ishandle(w)
    close(w); %waitbar
end
 
status = [status 'Tracking has finished.'];
setappdata(handles.main_figure,'status',status);



function yTopCrop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to yTopCrop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yTopCrop_edit as text
%        str2double(get(hObject,'String')) returns contents of yTopCrop_edit as a double
handles.yTopCrop = round(str2double(get(hObject,'String')));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yTopCrop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yTopCrop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function yBottomCrop_edit_Callback(hObject, eventdata, handles)
% hObject    handle to yBottomCrop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of yBottomCrop_edit as text
%        str2double(get(hObject,'String')) returns contents of yBottomCrop_edit as a double
handles.yBottomCrop = round(str2double(get(hObject,'String')));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function yBottomCrop_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to yBottomCrop_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in isCrop_checkbox.
function isCrop_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isCrop_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isCrop_checkbox
handles.isCrop = get(hObject,'Value');
guidata(hObject, handles);
