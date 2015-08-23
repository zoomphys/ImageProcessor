function outVars = getFRETZ(inVars)
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
%   startImg = This is the first image that we want to read in.  We will
%       begin to read in images starting with whichever imageis the
%       startImgth image in the folder.  This value will be an integer.
%       Spelled STARTIMG with only the I in uppercase and the rest of the
%       letters in lowercase.
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
if isfield(inVars,'isPassImages') && inVars.isPassImages
    Iy  = inVars.imYFP;
    Ic  = inVars.imCFP;
    startImg = 1;
else
    imdir = inVars.imdir;
    y_name = inVars.y_name;
    c_name = inVars.c_name;
    startImg = inVars.startImg;

    % Reads in the list of images we will use.
    yfiles = dir([imdir y_name]);
    cfiles = dir([imdir c_name]);
    Iy  = imread([imdir yfiles(startImg).name]);
    Ic  = imread([imdir cfiles(startImg).name]);
end

if isfield(inVars,'isImageParams') && inVars.isImageParams
    minArea = inVars.minArea;
    MATcrop = inVars.MATcrop;
    power   = inVars.power;
    numBits = inVars.numBits;
else
    % minArea = 50;
    minArea = 10;
    MATcrop = 50; %15;
    power   = 1.25;
    % croprad = 5;
    % emarg   = 1.5; %0.25;
    % bwcrppr = 20; %10;
    % frameps = 1;
    % cmprs   = 'none';
    numBits = 16;
    rDilate = 10;
    rErode = 1;
end

emptycell = struct;
masks{1} = [];

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
%--------------------------------------------------------------------------


%---------FINDING INITIAL CELL INFORMATION---------------------------------
% Reading in the images.
im  = Iy;

% keep track of raw images
Iyr = Iy;
Icr = Ic;

% Defines size of the images to come.
[r, c] = size(im);

% Finding the binary image.
im = im - min2(double(im));
im = double(im)/max2(double(im));
im = uint16(floor((2^numBits-1)*im));

% threshold image and get mask
bw = MovingAverageThresh(im, MATcrop, minArea, power);
tempim = im;
tempim(find(bw==0))=0;
bw = MovingAverageThresh(tempim, floor(MATcrop/2), minArea, power);

% Gets rid of the background noisefor Iy and Ic.
bw        = logical(bw);
backavgIy = mean(double(Iy(~bw)));
Iy        = uint16(floor(double(Iy) - backavgIy*ones(r, c)));
backavgIc = mean(double(Ic(~bw)));
Ic        = uint16(floor(double(Ic) - backavgIc*ones(r, c)));

% Finding the information of the initial cells.
L       = bwlabel(bw, 4);
stats   = regionprops(L, 'Area', 'Centroid');
stats12 = regionprops(imfill(L, 'holes'), 'Area', 'Centroid');
for i = 1:length(stats)
    bw = (L == i);
%     bw = logical(bw);
    bwt = logical(bw);
    bw = imerode(bwt, strel('disk', rErode, 0));
%     bw = imdilate(bwt, strel('disk', rDilate, 0));

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
outVars.tracks = tracks;
outVars.masks = masks;
outVars.L = imfill(L, 'holes');

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

