function output = gridFRET_single(inVars)
% Adapter from Darvin's FRETTracker
% Assume images are 16 bits

imdir = inVars.imdir;
y_name = inVars.y_name;
c_name = inVars.c_name;
startImg = inVars.startImg;

numRowDiv = inVars.numRowDiv;
numColDiv = inVars.numColDiv;
isThresholdFRET = inVars.isThresholdFRET;
isThresholdAdaptive = inVars.isThresholdAdaptive;
threshold_name = inVars.threshold_name;
threshold = inVars.threshold;

% Reads in the list of images we will use.
yfiles = dir(fullfile(imdir,y_name));
cfiles = dir(fullfile(imdir,c_name));
threshold_files = dir(fullfile(imdir,threshold_name));

% A lot of the hard-coded variables.
% for adaptive thresholding
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

%--------------------------------------------------------------------------


%---------FINDING INITIAL CELL INFORMATION---------------------------------
% In order to analyze the rest of the images and make a nice for-loop, we
% need to get the initial properties of the first frame.  We also start the
% troubleshooting movie here.

% Reading in the images.
Iy  = imread(fullfile(imdir,yfiles(startImg).name));
Ic  = imread(fullfile(imdir,cfiles(startImg).name));
Ithreshold  = imread(fullfile(imdir,threshold_files(startImg).name));

% keep track of the original raw images
Iyr = Iy;
Icr = Ic;

% Finding the binary image using YFP.
% im = Iy;

% Finding the binary image.
im = Ithreshold;
im = im - min2(double(im));
im = double(im)/max2(double(im));
im = uint16(floor((2^16-1)*im));

% Defines size of the images to come.
[r, c] = size(im);

% threshold image and get mask
if ~isThresholdFRET
    bw = ones(r,c);
else
    if isThresholdAdaptive
        % Finding the binary image.
        bw = MovingAverageThresh(im, MATcrop, minArea, power);
        % get part of image where the mask is 1
        tempim = im;
        tempim(find(bw==0))=0;
        bw = MovingAverageThresh(tempim, floor(MATcrop/2), minArea, power);
    else
        if threshold>0
            bw = im2bw(im,threshold);
        else
            bw = im2bw(im, graythresh(im));
        end
        bw = bwareaopen(bw, 4, 4);
    end
end

% Number of pixel in a division
pixRowDiv = floor(r/numRowDiv);
pixColDiv = floor(c/numColDiv);

bw        = logical(bw);

% Gets rid of the background noisefor Iy and Ic.
if ~isThresholdFRET
    backavgIy = 0;
    backavgIc = 0;
else
    backavgIy = mean(double(Iy(~bw)));
    backavgIc = mean(double(Ic(~bw)));
end

Iy        = uint16(floor(double(Iy) - backavgIy));
Ic        = uint16(floor(double(Ic) - backavgIc));

% Produce and image of the ratio c/y
Iratio = uint16(double(Ic)./double(Iy)*32768);
Iratio(~bw) = 0;

% Recording the total mask of the cells in question.
output.yfp = mean2(double(Iy(bw)));
output.cfp = mean2(double(Ic(bw)));
output.ratio = output.cfp/output.yfp;

output.rawy = mean2(double(Iyr(bw)));
output.rawc = mean2(double(Icr(bw)));
output.rawr = output.rawc/output.rawy;

output.yfpDiv = avgDiv(Iy,bw,numRowDiv,numColDiv);
output.cfpDiv = avgDiv(Ic,bw,numRowDiv,numColDiv);
output.ratioDiv = output.cfpDiv./output.yfpDiv;

output.rawyDiv = avgDiv(Iyr,bw,numRowDiv,numColDiv);
output.rawcDiv = avgDiv(Icr,bw,numRowDiv,numColDiv);
output.rawrDiv = output.rawcDiv./output.rawyDiv;

for i = 1:numRowDiv
    for j = 1:numColDiv
        iFlat = (i-1)*numColDiv+j;
        tracks(iFlat).x     = j;
        tracks(iFlat).y     = i;
        tracks(iFlat).frame = startImg;
        tracks(iFlat).area  = pixRowDiv*pixColDiv;
        tracks(iFlat).round = 0;
        tracks(iFlat).yfp   = output.yfpDiv(iFlat);
        tracks(iFlat).cfp   = output.cfpDiv(iFlat);
        tracks(iFlat).ratio = output.ratioDiv(iFlat);
        tracks(iFlat).rawy  = output.rawyDiv(iFlat);
        tracks(iFlat).rawc  = output.rawcDiv(iFlat);
        tracks(iFlat).rawr  = output.rawrDiv(iFlat);
    end
end


%--------------------------------------------------------------------------
% variables to output
output.tracks = tracks;
output.bw = bw;
output.Iratio = Iratio;

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
    function [output] = min2(input)
        output = min(min(input));
    end
%--------------------------------------------------------------------------


%---------FUNCTION: max2---------------------------------------------------
    function [output] = max2(input)
        output = max(max(input));
    end
%--------------------------------------------------------------------------


%---------FUNCTION: sum2---------------------------------------------------
    function [output] = sum2(input)
        output = sum(sum(input));
    end
%--------------------------------------------------------------------------


%---------FUNCTION: avgDiv---------------------------------------------------
    function [avgs] = avgDiv(im,mask,numRowDiv,numColDiv)
        [r, c] = size(im);
        pixRowDiv = floor(r/numRowDiv);
        pixColDiv = floor(c/numColDiv);

        avgs = [];
        for iRow=1:numRowDiv
            for iCol=1:numColDiv
                % make a mask of the division
                maskDiv = zeros(r,c);
                maskDiv((iRow-1)*pixRowDiv+1:iRow*pixRowDiv,(iCol-1)*pixColDiv+1:iCol*pixColDiv) = 1;
                cropped = double(im(mask&maskDiv));
                avgs = [avgs mean(cropped)];
            end
        end
    end
%--------------------------------------------------------------------------


end

