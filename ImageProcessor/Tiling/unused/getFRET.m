function output = getFRET(input)
% 1/20/2013: Adapted from Darwin's FRETtracker.m

%---------DESCRIPTION------------------------------------------------------
%FRETtracker Takes in a movie of FRET images of cells.  It then proceeds to
%track each cell, find its positions, velocities, areas, roundness, CFP
%intensities, YFP intensities, and the appropriate FRET ratios (CFP/YFP).
%   imdir    = This is the folder path that contains the images.  It is
%       assumed that the YFP and CFP intensity images are all in the same
%       folder.  This value will be a string.  Spelled IMDIR in all
%       lowercase.
%   y_name   = This is the name of all the YFP images.  It is recommended
%       that you use regex in this variable name.  This value will be a
%       string.  Spelled Y[UNDERSCORE]NAME in all lowercase.
%   c_name   = This is the name of all the CFP images.  It is recommended
%       that you use regex in this variable name.  This value will be a
%       string.  Spelled C[UNDERSCORE]NAME in all lowercase.
%   dt       = This is the time between each frame in minutes.  This value
%       will be used to find the exact velocities in terms of minutes.
%       This value will be a double.  Spelled DT in all lowercase.
%   startImg = This is the first image that we want to read in.  We will
%       begin to read in images starting with whichever imageis the
%       startImgth image in the folder.  This value will be an integer.
%       Spelled STARTIMG with only the I in uppercase and the rest of the
%       letters in lowercase.
%   numImg   = This is the number of images we will read in.  Starting with
%       the startImg, we will only read in numImg number of images.  This
%       value will be an integer.  Spelled NUMIMG with only the I in
%       capital and the rest of the letters in lowercase.
%   tau      = This is the amount of time that we will allow a cell to
%       disappear for before we give up on it.  For optimal results, you
%       should make tau equal to dt.  This value will be a double.  Spelled
%       TAU in all lowercase.
%   tau2     = This is basically the length of time we will allow a cell to
%       be still looked for "ish."  Note, it is very important to
%       differentiate this tau2 and tau (see above).  The variable tau2
%       only serves in the global assay of recombining tracks, not the
%       stalking.  Spelled TAU2 in all lowercase.
%   tracks   = This is basically the final goal of the program.  It will
%       contain all the information of each cell's tracks.  It will contain
%       the centroid of the cell at each frame it was found (x and y), the 
%       velocity of the cell at each frame it was found (vx and vy), the
%       area of the cell at each frame (area), the roundness of the cell at
%       each frame (round), the frames that the cell was found in (frame),
%       the background subtracted values of the YFP, CFP, and the CFP/YFP
%       ratio (yfp, cfp, and ratio respectively), and the raw values of the
%       YFP, CFP, and the CFP/YFP ratio (rawy, fawc, and rawr
%       respectively).  Spelled TRACKS in all lowercase.
%--------------------------------------------------------------------------


%---------SETTING VARIABLES------------------------------------------------
% This section will make some changes to predefined variables and also make
% other new variables from those predefined variables.  A lot of the time,
% it is just setting up some arrays or something along those lines.  Some
% things just should not be in the input arguments, as if we put everything
% into the input arguments, we'd be inputting input arguments all day long.
% Thus, the few hard-coded things will be here in the code.
% Fix the image directory name.

imdir = input.imdir;
y_name = input.y_name;
c_name = input.c_name;
startImg = input.startImg;

% Reads in the list of images we will use.
yfiles = dir([imdir y_name]);
cfiles = dir([imdir c_name]);

% A lot of the hard-coded variables.
minArea = 50;
emarg   = 1.5; %0.25;
MATcrop = 50; %15;
croprad = 5;
power   = 1.25;
bwcrppr = 20; %10;
frameps = 1;
cmprs   = 'none';

% This section will initialize the final variable "tracks."
emptycell.x     = [];
emptycell.y     = [];
emptycell.vx    = [];
emptycell.vy    = [];
emptycell.area  = [];
emptycell.round = [];
emptycell.frame = [];
emptycell.yfp   = [];
emptycell.cfp   = [];
emptycell.ratio = [];
emptycell.rawy  = [];
emptycell.rawc  = [];
emptycell.rawr  = [];
tracks(1)       = deal(emptycell);
totalcell.yfp   = [];
totalcell.cfp   = [];
totalcell.ratio = [];
totalcell.rawy  = [];
totalcell.rawc  = [];
totalcell.rawr  = [];
tottracks(1)    = deal(totalcell);
masks{1}        = [];
%--------------------------------------------------------------------------


%---------FINDING INITIAL CELL INFORMATION---------------------------------
% In order to analyze the rest of the images and make a nice for-loop, we
% need to get the initial properties of the first frame.  We also start the
% troubleshooting movie here.

% Reading in the images.
im  = imread([imdir yfiles(startImg).name]);
Iy  = im;
Ic  = imread([imdir cfiles(startImg).name]);
Iyr = Iy;
Icr = Ic;

% Defines size of the images to come.
[r, c] = size(im);

% Finding the binary image.
im = im - min2(double(im));
im = double(im)/max2(double(im));
im = uint16(floor((2^16-1)*im));
bw = MovingAverageThresh(im, MATcrop, minArea, power);

tempim = im;
for i = 1:r
    for j = 1:c
        if bw(i, j) == 0
            tempim(i, j) = 0;
        end
    end
end
bw = MovingAverageThresh(tempim, floor(MATcrop/2), minArea, power);

% Gets rid of the background noisefor Iy and Ic.
bw        = logical(bw);
backavgIy = mean(double(Iy(~bw)));
Iy        = uint16(floor(double(Iy) - backavgIy*ones(r, c)));
backavgIc = mean(double(Ic(~bw)));
Ic        = uint16(floor(double(Ic) - backavgIc*ones(r, c)));

% Recording the total mask of the cells in question.
tottracks.yfp   = mean2(double(Iy(bw)));
tottracks.cfp   = mean2(double(Ic(bw)));
tottracks.ratio = tottracks.cfp/tottracks.yfp;
tottracks.rawy  = mean2(double(Iyr(bw)));
tottracks.rawc  = mean2(double(Icr(bw)));
tottracks.rawr  = tottracks.rawc/tottracks.rawy;

% Finding the information of the initial cells.
L       = bwlabel(bw, 4);
stats   = regionprops(L, 'Area', 'Centroid');
stats12 = regionprops(imfill(L, 'holes'), 'Area', 'Centroid');
for i = 1:length(stats)
    bw = (L == i);
    bw = logical(bw);
    tracks(i).x     = stats(i).Centroid(1);
    tracks(i).y     = stats(i).Centroid(2);
    tracks(i).vx    = 0;
    tracks(i).vy    = 0;
    tracks(i).frame = startImg;
%     tracks(i).frame = 1;
    tracks(i).area  = stats12(i).Area;
    tracks(i).round = 4*pi*stats12(i).Area/sum2(double(bwperim(bw, 8)));
    tracks(i).yfp   = mean(double(Iy(bw)));
    tracks(i).cfp   = mean(double(Ic(bw)));
    tracks(i).ratio = tracks(i).cfp/tracks(i).yfp;
    tracks(i).rawy  = mean(double(Iyr(bw)));
    tracks(i).rawc  = mean(double(Icr(bw)));
    tracks(i).rawr  = tracks(i).rawc/tracks(i).rawy;
    masks{i} = bw;
end


%--------------------------------------------------------------------------
% variables to output
output.tracks = tracks;
output.tottracks = tottracks;
output.masks = masks;
output.L = imfill(L, 'holes');

%---------FUNCTION: MovingAverageThresh------------------------------------
    function [binaryImg] = MovingAverageThresh(Img, Cropper, MINarea, Power)
        % Defining the average image filter.
        hfilter = fspecial('average', Cropper*2 + 1);
        
        % Finding the doubly average subtracted image.
        tempImg = imfilter(Img, hfilter, 'replicate');
        tempImg = imfilter(tempImg, hfilter, 'replicate');
        tempImg = Img - tempImg;
        
        % Raising the image to a power to exagerate the difference.
        tempImg = uint16(floor(double(tempImg).^Power));
        
        % Making the image binary.
        binaryImg = im2bw(tempImg, graythresh(tempImg));
        binaryImg = bwareaopen(binaryImg, MINarea, 4);
    end        
%--------------------------------------------------------------------------


%---------FUNCTION: min2---------------------------------------------------
    function [matrixMIN] = min2(inputMIN)
        matrixMIN = min(min(inputMIN));
    end
%--------------------------------------------------------------------------


%---------FUNCTION: max2---------------------------------------------------
    function [matrixMAX] = max2(inputMAX)
        matrixMAX = max(max(inputMAX));
    end
%--------------------------------------------------------------------------


%---------FUNCTION: sum2---------------------------------------------------
    function [matrixSUM] = sum2(inputSUM)
        matrixSUM = sum(sum(inputSUM));
    end
%--------------------------------------------------------------------------


%---------GLOSSERY---------------------------------------------------------
% a
% aviobj
% backavgIc
% backavgIy
% bw
% bw_crop
% bwt
% bwcrppr
% c
% c1
% c1lower
% c1upper
% cfiles
% circ
% cmap
% cmprs
% comp
% comp1
% counter
% croprad
% cropxh
% cropxl
% cropyh
% cropyl
% cx
% cy
% dist
% emarg
% emptycell
% frameps
% Ic
% Icr
% Iy
% Iyr
% i
% i1
% im
% im_crop
% j
% j1
% k
% L
% MATcrop
% masks
% maxf
% maxf2
% maxl
% maxt
% message
% minArea
% modder
% power
% r
% r1
% r1lower
% r1upper
% rad
% stats
% stats12
% t
% tempim
% temptracks
% totalmasks
% w
% yfiles
% MovingAverageThresh
%   binaryImg
%   Cropper
%   hfilter
%   Img
%   MINarea
%   Power
%   tempImg
% max2
%   inputMAX
%   matrixMAX
% min2
%   inputMIN
%   matrixMIN
% sum2
%   inputSUM
%   matrixSUM
%--------------------------------------------------------------------------


%---------HISTORY----------------------------------------------------------
% Created by Darvin Yi.                                          2012-08-07
% Modified by Darvin Yi.                                         2012-08-07
%--------------------------------------------------------------------------


%---------INDEX------------------------------------------------------------
% DESCRIPTION
% SETTING VARIABLES
% FINDING INITIAL CELL INFORMATION
% ANALYZING ALL OTHER FRAMES TO CONNECT CELLS
% FUNCTION: MovingAverageThresh
% GLOSSERY
% HISTORY
% INDEX
%--------------------------------------------------------------------------
end

