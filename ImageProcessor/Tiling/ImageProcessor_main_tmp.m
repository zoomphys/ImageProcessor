% Author: Duong "Zoom" Nguyen
% Email: duongnh "at" gmail "dot" com
% version 1.7
% Jan 22, 2015

function varargout = ImageProcessor_main_tmp(varargin)
% IMAGEPROCESSOR_MAIN_TMP MATLAB code for ImageProcessor_main_tmp.fig
%      IMAGEPROCESSOR_MAIN_TMP, by itself, creates a new IMAGEPROCESSOR_MAIN_TMP or raises the existing
%      singleton*.
%
%      H = IMAGEPROCESSOR_MAIN_TMP returns the handle to a new IMAGEPROCESSOR_MAIN_TMP or the handle to
%      the existing singleton*.
%
%      IMAGEPROCESSOR_MAIN_TMP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in IMAGEPROCESSOR_MAIN_TMP.M with the given input arguments.
%
%      IMAGEPROCESSOR_MAIN_TMP('Property','Value',...) creates a new IMAGEPROCESSOR_MAIN_TMP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ImageProcessor_main_tmp_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ImageProcessor_main_tmp_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ImageProcessor_main_tmp

% Last Modified by GUIDE v2.5 28-Jan-2015 08:57:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ImageProcessor_main_tmp_OpeningFcn, ...
                   'gui_OutputFcn',  @ImageProcessor_main_tmp_OutputFcn, ...
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


% --- Executes just before ImageProcessor_main_tmp is made visible.
function ImageProcessor_main_tmp_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ImageProcessor_main_tmp (see VARARGIN)

libdir = fullfile(pwd,'lib');
if exist(libdir,'dir')
    addpath(genpath(libdir));
end

% Choose default command line output for ImageProcessor_main_tmp
handles.output = hObject;

%%%%% Initialize variables
handles.filedir = get(handles.filedir_edit,'String');
handles.prefix = get(handles.prefix_edit,'String');
handles.dirQueue = {};

handles.isAutoContrast = get(handles.isAutoContrast_checkbox,'Value');

handles.isLog = get(handles.isLog_checkbox,'Value');

handles.isShowRightImage = 1;
handles.maxPixelValue = 65535;

% Variables to save to a file so that they are loaded at the next startup
handles.mainVars = struct;

handles.mainVars.guivarsFileName = get(handles.varsFileName_edit,'String');
handles.mainVars.savedVars = {
    'filedir' 'prefix' 'suffix' 'imgExtension' 'prefixQueue' ...
    'isAutoContrast' 'isLog' 'isShowRightImage' ...
    };
handles.mainVars.guiName = {
    'filedir_edit' 'prefix_edit' 'suffix_edit' 'imgExtension_edit' 'prefixQueue_edit' ...
    'isAutoContrast_checkbox' 'isLog_checkbox' 'isShowRightImage_checkbox' ...
    };
handles.mainVars.guiType = {
    'String' 'String' 'String' 'String' 'String' ...
    'Value' 'Value' 'Value' ...
    };
handles = getVars(handles,handles.mainVars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize status and log
setappdata(handles.main_figure,'log',cell(0));
setappdata(handles.main_figure,'status','');

% Make functions public
handles.loadFiles = @loadFiles;
handles.load_pushbutton_Callback = @load_pushbutton_Callback;
handles.current_edit_Callback = @current_edit_Callback;

% Update handles structure
updatestatus(handles,'Welcome!');
guidata(hObject, handles);

% UIWAIT makes ImageProcessor_main_tmp wait for user response (see UIRESUME)
% uiwait(handles.main_figure);


% --- Outputs from this function are returned to the command line.
function varargout = ImageProcessor_main_tmp_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in exit_pushbutton.
function exit_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to exit_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
main_figure_CloseRequestFcn(handles.main_figure, eventdata, handles);
% delete(hObject);


% --- Executes when user attempts to close main_figure.
function main_figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to main_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
promptMessage = sprintf('Do you want to Continue exiting?\n(Click Cancel to stay running)');
selectedButton = questdlg(promptMessage, 'Exit Dialog','Continue exiting', 'Cancel', 'Continue exiting');
if strcmp(selectedButton, 'Cancel')        % Stay in the program. Do not exit.
    return;
end

% save gui variables
handles = saveVars(handles,handles.mainVars);
if isfield(handles,'isSubFigureOpened') && handles.isSubFigureOpened
    handles = saveVars(handles,handles.subVars);
end    

updatestatus(handles,'Goodbye!');

% save log
if handles.isLog
    saveLog_pushbutton_Callback(hObject, eventdata, handles);
end

libdir = fullfile(pwd,'lib');
if exist(libdir,'dir')
    rmpath(genpath(libdir));
end

% Continue to exit by deleting this GUI
delete(hObject);


% --- Executes during object creation, after setting all properties.
function filedir_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filedir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in browse_pushbutton.
function browse_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to browse_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dirstr = uigetdir(handles.filedir,'Specify the first image folder...');
if ~ischar(dirstr)
    updatestatus(handles,'Not a valid path.');
    return;
end

%Always use / to separate folders
%No ending slash
handles.filedir = [strrep(dirstr,'\','/') '/'];
set(handles.filedir_edit,'String',handles.filedir);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function current_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function prefix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%%%%%%%%%%%%%%%%%%%%%%% Begin user-defined functions

%%%%%%%%%%%%%%%%%%%%%%% End user-defined functions


function handles = loadImageInSeries(num,handles)
% load an image in a series specified by num

% check if num is integer and in range
if mod(num,1)~=0 || num<1 || num>handles.numFiles
    updatestatus(handles,'Specified index is out of range.');
    return;
end

% set(handles.fileName_text,'String',handles.fileNames{num});

currentFilePath = fullfile(handles.filedir, handles.fileNames{num});
handles = loadCurrentImage(handles,currentFilePath);

if ~isempty(handles.status)
     updatestatus(handles,handles.status);
     return;
end


function handles = current_edit_Callback(hObject, eventdata, handles)
% hObject    handle to current_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of current_edit as text
%        str2double(get(hObject,'String')) returns contents of current_edit as a double
num = round(str2double(get(hObject,'String')));
handles = loadImageInSeries(num,handles);
handles.iCurrentFile = num;

set(handles.current_slider,'Value',num);

guidata(hObject, handles);


function numfiles_edit_Callback(hObject, eventdata, handles)
% hObject    handle to numfiles_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numfiles_edit as text
%        str2double(get(hObject,'String')) returns contents of numfiles_edit as a double
num = str2double(get(hObject,'String'));
% check if num is integer and in range
if mod(num,1)~=0 || num<1 || num>handles.numFiles-handles.iCurrentFile+1
    updatestatus(handles,'Specified number if images is not valid.');
    return;
end
handles.numFiles = num;
updatestatus(handles,'Number of images has been updated.');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function numfiles_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numfiles_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function filedir_edit_Callback(hObject, eventdata, handles)
% hObject    handle to filedir_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filedir_edit as text
%        str2double(get(hObject,'String')) returns contents of filedir_edit as a double
filedir = get(handles.filedir_edit,'String');

status = 1;
if ~exist(filedir,'dir')
    promptMessage = sprintf('Directory does not exist. Do you want to create it?\n(Click Cancel to stay running)');
    selectedButton = questdlg(promptMessage, 'Create Directory Dialog','Create', 'Cancel', 'Cancel');
    if strcmp(selectedButton, 'Create')
        status = mkdir(filedir);
        if ~status
            updatestatus(handles,'Cannot create directory.');
        end
    else
        status = 0;
    end
end

if ~status
    set(handles.filedir_edit,'String',handles.filedir);
    updatestatus(handles,'Invalid directory.');
    return;
end

filedir = strrep(filedir,'\','/');
if filedir(end) ~= '/'
    handles.filedir = [filedir '/'];
else
    handles.filedir = filedir;
end

set(handles.filedir_edit,'String',handles.filedir);
updatestatus(handles,'File directory has been updated.');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function left_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to left_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate left_axes
set(0,'defaultaxeslinewidth',1);


% --- Executes during object creation, after setting all properties.
function right_axes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to right_axes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate right_axes
set(0,'defaultaxeslinewidth',1);


function handles = loadCurrentImage(handles,currentFilePath)
% Read and display image from currentFilePath

handles.status = '';

try
   imgTemp = imread(currentFilePath);
catch exception
    handles.status = ['Error reading image from file' currentFilePath];
    return;
end

if ~strcmp(handles.colorType,'grayscale')
   imgTemp = uint16(rgb2gray(imgTemp));
   handles.bitDepth = 16;
end

if handles.isAutoContrast
    handles.imgLeft = autocontrast(imgTemp,handles.bitDepth);
else
    handles.imgLeft = imgTemp;
end
    
%axes(handles.left_axes);
dispImage(handles.imgLeft,handles);


% --- Executes on button press in load_pushbutton.
function handles = load_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to load_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filePattern = [handles.prefix '*' handles.suffix '.' handles.imgExtension];
filePattern 
outVars = loadFiles(handles,handles.filedir,filePattern);
if ~isstruct(outVars)
     updatestatus(handles,'Error loading files.');
     return;
end

handles.iCurrentFile = 1;
handles.files = outVars.files;
handles.fileNames = outVars.fileNames;
handles.numFiles = outVars.numFiles;
handles.currentFilePath = outVars.firstFilePath;
handles.fileInfo = outVars.fileInfo;
handles.bitDepth = outVars.bitDepth;
handles.colorType = outVars.colorType;

handles = loadCurrentImage(handles,handles.currentFilePath);
if ~isempty(handles.status)
     updatestatus(handles,handles.status);
     return;
end

% Initialize detected cell information
if isfield(handles,'numSpots')
    handles = rmfield(handles,'numSpots');
end

set(handles.current_edit,'String',num2str(handles.iCurrentFile));
set(handles.numfiles_text,'String',['/  ',num2str(handles.numFiles)]);

if handles.numFiles==1
    maxSlider = 2; % so that the max value is bigger than the min value
else
    maxSlider = handles.numFiles;
end
set(handles.current_slider,'Value',1,'Max',maxSlider,'SliderStep',[1/(maxSlider-1) 1/(maxSlider-1)]);

updatestatus(handles,['Success loading files from: ' handles.filedir]);

guidata(hObject, handles);


function outVars = loadFiles(handles,filedir,filePattern)
% load all images given by filePattern in directory filedir
files = dir(fullfile(filedir,filePattern));
fileNames = {files(:).name};
numFiles = length(fileNames);

if numFiles<1
    updatestatus(handles,'No file is found.');
    outVars = nan;
    return;
end

firstFilePath = fullfile(filedir, fileNames{1});
try
   fileInfo = imfinfo(firstFilePath);
catch exception
    updatestatus(handles,['Error reading image information from file' firstFilePath]);
    return;
end

if isfield(fileInfo,'BitDepth')
    bitDepth = uint8(fileInfo.BitDepth);
else
    bitDepth = 16;
end
if isfield(fileInfo,'ColorType')
    colorType = fileInfo.ColorType;
else
    colorType = 'grayscale';
end

outVars.files = files;
outVars.fileNames = fileNames;
outVars.numFiles = numFiles;
outVars.firstFilePath = firstFilePath;
outVars.fileInfo = fileInfo;
outVars.bitDepth = bitDepth;
outVars.colorType = colorType;


function prefix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to prefix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefix_edit as text
%        str2double(get(hObject,'String')) returns contents of prefix_edit as a double
handles.prefix = get(hObject,'String');
guidata(hObject, handles);


% --- Executes on slider movement.
function current_slider_Callback(hObject, eventdata, handles)
% hObject    handle to current_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
num = get(hObject,'Value');
if handles.numFiles==1
    num_index = 1;
else
    num_index = round(num);
end

handles = loadImageInSeries(num_index,handles);
handles.iCurrentFile = num_index;

set(handles.current_edit,'String',num2str(num_index));
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function current_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in isAutoContrast_checkbox.
function isAutoContrast_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isAutoContrast_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isAutoContrast_checkbox
handles.isAutoContrast = get(hObject,'Value');
guidata(hObject, handles);


function varsFileName_edit_Callback(hObject, eventdata, handles)
% hObject    handle to varsFileName_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of varsFileName_edit as text
%        str2double(get(hObject,'String')) returns contents of varsFileName_edit as a double
handles.varsFileName = get(hObject,'String');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function varsFileName_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to varsFileName_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in isLog_checkbox.
function isLog_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isLog_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isLog_checkbox
handles.isLog = get(hObject,'Value');
guidata(hObject, handles);




% --- Executes on button press in isShowRightImage_checkbox.
function isShowRightImage_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to isShowRightImage_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of isShowRightImage_checkbox
handles.isShowRightImage = get(hObject,'Value');
guidata(hObject, handles);


function suffix_edit_Callback(hObject, eventdata, handles)
% hObject    handle to suffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of suffix_edit as text
%        str2double(get(hObject,'String')) returns contents of suffix_edit as a double
handles.suffix = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function suffix_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to suffix_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imgExtension_edit_Callback(hObject, eventdata, handles)
% hObject    handle to imgExtension_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imgExtension_edit as text
%        str2double(get(hObject,'String')) returns contents of imgExtension_edit as a double
handles.imgExtension = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function imgExtension_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imgExtension_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function dispImage(img,handles)
% display image on the axes
if handles.isAutoContrast
    img = autocontrast(img,handles.bitDepth);
end 
imshow(img,'Parent',handles.left_axes);


% --- Executes on button press in saveLeftImg_pushbutton.
function saveLeftImg_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveLeftImg_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savedImgPath =  [handles.filedir 'imgLeft_' getTime_filename() '.png'];
%img = getframe(handles.left_axes);
%imwrite(img.cdata,savedImgPath);
export_fig(savedImgPath,handles.left_axes);
updatestatus(handles,['Left image has been saved to file: ' savedImgPath]);


% --- Executes on button press in saveRightImg_pushbutton.
function saveRightImg_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveRightImg_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savedImgPath =  [handles.filedir 'imgRight_' getTime_filename() '.png'];
% img = getframe(handles.right_axes);
% imwrite(img.cdata,savedImgPath);
export_fig(savedImgPath,handles.right_axes);
updatestatus(handles,['Right image has been saved to file: ' savedImgPath]);


% --- Executes on button press in saveLog_pushbutton.
function saveLog_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveLog_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~exist(handles.filedir,'dir')
    filedir = '';
    updatestatus(handles,'Directory does not exist. Log is saved to the program directory.');
else
    filedir = handles.filedir;
end

log = getappdata(handles.main_figure,'log');
logFilePath = fullfile(filedir,'ImageProcessor.log');

updatestatus(handles,['Current log is being saved to: ' logFilePath]);
fid = fopen(logFilePath,'at');
fprintf(fid, '%s\n', log{:});
fclose(fid);
setappdata(handles.main_figure,'log',cell(0));


% --- Executes on button press in saveFigure_pushbutton.
function saveFigure_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveFigure_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
savedImgPath =  [handles.filedir 'img_' getTime_filename() '.png'];
export_fig(savedImgPath,handles.main_figure);
updatestatus(handles,['Snapshot of gui has been saved to file: ' savedImgPath]);



% --- Executes on button press in showMenu_pushbutton.
function showMenu_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to showMenu_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentState = get(handles.main_figure,'MenuBar');
if strcmp(currentState,'figure')
    set(handles.main_figure,'MenuBar','None');  
    set(handles.showMenu_pushbutton,'String','Show Menu');
    updatestatus(handles,'Figure menu is now hidden.');
else
    set(handles.main_figure,'MenuBar','figure');  
    set(handles.showMenu_pushbutton,'String','Hide Menu');
    updatestatus(handles,'Figure menu is now shown.');
end


% --- Executes on button press in dictyTracker_pushbutton.
function dictyTracker_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to dictyTracker_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.subGUIName = 'DictyTracker';
addpath(handles.subGUIName);

h = openfig(handles.subGUIName,'new','invisible');
handles.subGUIHandles = get(h,'children');
set(handles.subGUIHandles ,'Parent',handles.addon_uipanel);
delete(h);

handles.subGUITags = get(handles.subGUIHandles,'Tag');
handles.subGUINumTags = length(handles.subGUITags);
handles.subGUITransfered = logical(ones(1,handles.subGUINumTags));
for i=1:handles.subGUINumTags
    var = handles.subGUITags{i};
    if isfield(handles,var)
        updatestatus(handles,['SubGUI tag ',var,' already exists in the main GUI.']);
        handles.subGUITransfered(i) = false;
    else
        handles.(var) = handles.subGUIHandles(i);
    end
end


handles = init_guivars(handles,'DictyTracker_vars.xlsx');
handles.subVars = handles.outVars;
handles = getVars(handles,handles.subVars);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables specific to this sub figure
handles.savedFields = handles.subVars.savedFields;
set(handles.plotFields_popupmenu,'String',handles.savedFields);

% Color map for plotting multiple curves
handles.colors = lines(10);

handles.isSubFigureOpened = true;
updatestatus(handles,'DictyTracker SubGUI has been opened.');

guidata(hObject, handles);

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function main_figure_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to dicty_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mousepos = get(handles.left_axes,'Currentpoint');
xpos = mousepos(1,1)-1;
ypos = mousepos(1,2)-1;

if ~isfield(handles,'imgLeft')
    updatestatus(handles,'No image has been loaded to the left axes.');
    return;
end

[width,height] = size(handles.imgLeft);
if xpos<1 || xpos>width
    updatestatus(handles,'X position is out of range.');
    return;
end
if ypos<1 || ypos>height
    updatestatus(handles,'Y position is out of range.');
    return;
end
     
% convert y coordinate from image to axes. In matlab these two go different ways
%ypos = handles.imheight - ypos;

% % Now we find out which mouse button was clicked, and whether a
% % keyboard modifier was used, e.g. shift or ctrl
switch get(gcf,'SelectionType')
  case 'normal' % Click left mouse button.
        xpos = round(xpos);
        ypos = round(ypos);
        if isfield(handles,'isManualSelect') && handles.isManualSelect
            handles.selectedPosArray = [handles.selectedPosArray;[xpos,ypos]];
            updatestatus(handles,['Current position: (x,y,value) = (' num2str(xpos) ',' num2str(ypos) ',' num2str(handles.imgLeft(ypos,xpos)) ') has been added to the position array.']);
        else
            updatestatus(handles,['Current position: (x,y,value) = (' num2str(xpos) ',' num2str(ypos) ',' num2str(handles.imgLeft(ypos,xpos)) ').']);
        end
        
  case 'alt' % Ctrl - click left mouse button or click right mouse button.
        xpos = round(xpos);
        ypos = round(ypos);
        handles.selectedSpot = [xpos,ypos];
        updatestatus(handles,['Current position: (x,y,value) = (' num2str(xpos) ',' num2str(ypos) ',' num2str(handles.imgLeft(ypos,xpos)) ') has been added to the position array.']);
        
  case 'extend' % Shift - click left mouse button or click both left and right mouse buttons.
        % find nearest spot near mouse click and toggle the selection
        if ~isfield(handles,'allSpotPos')
            updatestatus(handles,'No cell has been detected.');
            return;
        end
        handles = selectSpot(handles,[xpos,ypos]);
        
  case 'open'   % Double-click any mouse button.
end
guidata(hObject, handles);


function handles = selectSpot(handles,position)
% Select nearest cell to position where position is [x,y]
index = nearest(handles.allSpotPos,position);
handles.selectedSpot = handles.allSpotPos(index,:);

if handles.selectedIndexArray(index)
    symbol = 'r+';
    newStatus = 'deselected';
else
    symbol = 'g+';
    newStatus = 'selected';
end

handles.selectedIndexArray(index) = 1-handles.selectedIndexArray(index);
hold(handles.left_axes,'on');
plot(handles.left_axes,handles.selectedSpot(1),handles.selectedSpot(2),symbol);
hold(handles.left_axes,'off');

posStr = ['(' num2str(round(handles.selectedSpot(1))) ',' num2str(round(handles.selectedSpot(2))) ')'];
updatestatus(handles,['Spot number ' num2str(index) ' at position ' posStr ' has been ' newStatus '.']);


% --- Executes on button press in manualSelect_togglebutton.
function manualSelect_togglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to manualSelect_togglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of manualSelect_togglebutton
handles.isManualSelect = get(hObject,'Value');
if handles.isManualSelect
    handles.selectedPosArray = [];
    updatestatus(handles,'Manual selection is now on. Please click on the objects on the left image.');
else
    updatestatus(handles,'Manual selection is now off.');
end

guidata(hObject, handles);


% --- Executes when main_figure is resized.
function main_figure_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to main_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in closeSubGUI_pushbutton.
function closeSubGUI_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to closeSubGUI_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'isSubFigureOpened') || ~handles.isSubFigureOpened
    updatestatus(handles,'No SubGUI is opened.');
    return;
end

handles = saveVars(handles,handles.subVars);

for i=1:handles.subGUINumTags
    delete(handles.subGUIHandles(i));
    
    var = handles.subGUITags{i};
    if handles.subGUITransfered(i)
        if isfield(handles,var)
            handles = rmfield(handles,var);
%             disp(var)
        else
            updatestatus(handles,['SubGUI tag ',var,' does not exist in the main GUI.']);
        end
    end
end
        
rmpath(handles.subGUIName);
handles = rmfield(handles,{'subGUIHandles','subGUITags','subGUINumTags','subGUITransfered','subGUIName'});

updatestatus(handles,'SubGUI has been closed.');
handles.isSubFigureOpened = false;

guidata(hObject, handles);


% --- Executes on button press in addToQueue_pushbutton.
function addToQueue_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addToQueue_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dirQueue{end+1} = handles.filedir;
updatestatus(handles,['Current queue: ' handles.dirQueue]);
guidata(hObject, handles);


function prefixQueue_edit_Callback(hObject, eventdata, handles)
% hObject    handle to prefixQueue_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of prefixQueue_edit as text
%        str2double(get(hObject,'String')) returns contents of prefixQueue_edit as a double
handles.prefixQueue = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function prefixQueue_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to prefixQueue_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in addQueueSeries_pushbutton.
function addQueueSeries_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to addQueueSeries_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pattern = fullfile(handles.filedir,[handles.prefixQueue '*']);
files = dir(pattern);
files = files(find([files.isdir]));

if length(files)<1
    updatestatus(handles,'No subdirs found.');
    return;
end
for iDir=1:length(files)
    handles.dirQueue{end+1} = fullfile(handles.filedir,files(iDir).name);
end

updatestatus(handles,['Current queue: ' handles.dirQueue]);
guidata(hObject, handles);


% --- Executes on button press in tiling_pushbutton.
function tiling_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to tiling_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.subGUIName = 'Tiling';
addpath(handles.subGUIName);

% open new figure for new gui and embed its elements into addon uipanel
h = openfig(handles.subGUIName,'new','invisible');
handles.subGUIHandles = get(h,'children');
set(handles.subGUIHandles ,'Parent',handles.addon_uipanel);
delete(h);

% record subgui info uicontrols and keep track of the names in conflict
handles.subGUITags = get(handles.subGUIHandles,'Tag');
handles.subGUINumTags = length(handles.subGUITags);
handles.subGUITransfered = logical(ones(1,handles.subGUINumTags));
for i=1:handles.subGUINumTags
    var = handles.subGUITags{i};
    if isfield(handles,var)
        updatestatus(handles,['SubGUI tag ',var,' already exists in the main GUI.']);
        handles.subGUITransfered(i) = false;
    else
        handles.(var) = handles.subGUIHandles(i);
    end
end

handles = init_guivars(handles,'Tiling_vars.xlsx');
handles.subVars = handles.outVars;

handles.profiles = handles.subVars.profiles;
set(handles.profiles_popupmenu,'String', handles.profiles);

handles = getVars(handles,handles.subVars);

handles.isSubFigureOpened = true;
updatestatus(handles,'Tiling SubGUI has been opened.');

guidata(hObject, handles);


% --- Executes on button press in clearQueue_pushbutton.
function clearQueue_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to clearQueue_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dirQueue = {};
updatestatus(handles,'Queue has been cleared');
guidata(hObject, handles);
