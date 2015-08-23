function [tracks, tottracks] = FRETtracker(rootdir, y_name, c_name, yy_name, dt, startImg, numImg, tau, tau2, positions, alpha, beta)
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
%----For each position in the experiment-----------------------------------
 for p = 1:positions;
    
    stringPosition = sprintf('Now working on position %d', p);
    disp(stringPosition);
    posPrefix = 'Pos';
    suffix = num2str(p-1);  
    subFolder = strcat(posPrefix,suffix);

%---------SETTING VARIABLES------------------------------------------------
% This section will make some changes to predefined variables and also make
% other new variables from those predefined variables.  A lot of the time,
% it is just setting up some arrays or something along those lines.  Some
% things just should not be in the input arguments, as if we put everything
% into the input arguments, we'd be inputting input arguments all day long.
% Thus, the few hard-coded things will be here in the code.
% Fix the image directory name.
% Fix the image directory name.
imdir = [rootdir '\' subFolder '\'];

% Reads in the list of images we will use.
yfiles = dir([imdir y_name]);
cfiles = dir([imdir c_name]);
yyfiles = dir([imdir yy_name]);

% Converts tau into frames.
maxf  = tau/dt;
maxf2 = tau2/dt;

% A lot of the hard-coded variables.
xtalk   = 0;
minArea = 50;
emarg   = .33;
MATcrop = 50;
croprad = 5;
power   = 1.25;
bwcrppr = 25;
frameps = 4;
cmprs   = 'none';
erode   = 3;

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
emptycell.yyfp  = [];
emptycell.ratio = [];
emptycell.pratio = [];
emptycell.rawy  = [];
emptycell.rawc  = [];
emptycell.rawyy = [];
emptycell.rawr  = [];
emptycell.prawr = [];
emptycell.corr = [];
emptycell.pcorr = [];
emptycell.rawcorr = [];
emptycell.prawcorr = [];
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
% Setting up the video.
aviobj = avifile ([imdir 'FRETmovie.avi'], 'compression', cmprs, 'fps', frameps);

% Setting up the waitbar.
w = waitbar(0, 'Please wait...');
message = ['Analyzing Frame ' num2str(1) ' of ' num2str(numImg) ':'];
waitbar((1/numImg), w, message);

% Reading in the images.
im  = imread([imdir yyfiles(startImg).name]);
Iyy = imread([imdir yyfiles(startImg).name]);
Iy  = imread([imdir yfiles(startImg).name]);
Ic  = imread([imdir cfiles(startImg).name]);
Iyyr = imread([imdir yyfiles(startImg).name]);
Iyr = imread([imdir yfiles(startImg).name]);
Icr = imread([imdir cfiles(startImg).name]);

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
backavgIyy = mean(double(Iyy(~bw)));
Iyy       = uint16(floor(double(Iyy) - backavgIyy*ones(r, c)));

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
    bwt = logical(bw);
    bw = imerode(bwt, strel('disk', erode, 0));
    tracks(i).x     = stats(i).Centroid(1);
    tracks(i).y     = stats(i).Centroid(2);
    tracks(i).vx    = 0;
    tracks(i).vy    = 0;
    tracks(i).frame = 1;
    tracks(i).area  = stats12(i).Area;
    tracks(i).round = 4*pi*stats12(i).Area/sum2(double(bwperim(bw, 8)));
    tracks(i).yfp   = mean(double(Iy(bw)));
    tracks(i).cfp   = mean(double(Ic(bw)));
    tracks(i).yyfp  = mean(double(Iyy(bw)));
    tracks(i).ratio = tracks(i).cfp/tracks(i).yfp;
    tracks(i).pratio = mean(double(Ic(bw))./double(Iy(bw)));
    tracks(i).rawy  = mean(double(Iyr(bw)));
    tracks(i).rawc  = mean(double(Icr(bw)));
    tracks(i).rawyy = mean(double(Iyyr(bw)));
    tracks(i).rawr  = tracks(i).rawc/tracks(i).rawy;
    tracks(i).prawr = mean(double(Icr(bw))./double(Iyr(bw)));
    tracks(i).corr = (tracks(i).yfp - alpha*tracks(i).yyfp - beta*tracks(i).cfp)/(alpha*tracks(i).yyfp);
    tracks(i).pcorr = mean((double(Iy(bw)) - alpha*double(Iyy(bw)) - beta*double(Ic(bw)))/(alpha*double(Iyy(bw))));
    tracks(i).rawcorr = (tracks(i).rawy - alpha*tracks(i).rawyy - beta*tracks(i).rawc)/(alpha*tracks(i).rawyy);
    tracks(i).prawcorr = mean((double(Iyr(bw)) - alpha*double(Iyyr(bw)) - beta*double(Icr(bw)))/(alpha*double(Iyyr(bw))));
    masks{i} = bw;
end

% Adding the first frame to the movie.
cmap = hsv(64); %hsv(length(tracks));
im = double(im)/max2(double(im));
im = repmat(im, [1 1 3]);
totalmasks = zeros(r, c);
totalmasks = repmat(totalmasks, [1 1 3]);
for i = 1:length(tracks)
    rem(i,63)+1
    for t = 1:3
        totalmasks(:, :, t) = totalmasks(:, :, t) + cmap(rem(i,63)+1, t)*masks{i};
    end
end
aviobj = addframe(aviobj, [im totalmasks]);
modder = 64; %length(tracks);
%--------------------------------------------------------------------------


%---------ANALYZING ALL OTHER FRAMES TO CONNECT CELLS----------------------
% We will now go through all the other images.  At each image, we'll read
% in the corresponding image, and then do analysis on it.  First, we'll
% stalk the previous cells we found.  Then, we'll add cells that we found
% new.
for i = 2:numImg
    % Displaying the waitbar.
    message = ['Analyzing Frame ' num2str(i) ' of ' num2str(numImg) ':'];
    waitbar((i/numImg), w, message);
    
    % Reading in the images.
    im  = imread([imdir yyfiles(i-1+startImg).name]);
    Iyy = imread([imdir yyfiles(i-1+startImg).name]);
    Iyyr = imread([imdir yyfiles(i-1+startImg).name]);
    Iy  = imread([imdir yfiles(i-1+startImg).name]);
    Ic  = imread([imdir cfiles(i-1+startImg).name]);
    Iyr = imread([imdir yfiles(i-1+startImg).name]);
    Icr = imread([imdir cfiles(i-1+startImg).name]);
    
    % Saving Memory
    clear masks
    masks{1} = [];
    
    % Finds the binary image.
    im = im - min2(double(im));
    im = double(im)/max2(double(im));
    im = uint16(floor((2^16-1)*im));
    bw = MovingAverageThresh(im, MATcrop, minArea, power);
    tempim = im;
    for i1 = 1:r
        for j1 = 1:c
            if bw(i1, j1) == 0
                tempim(i1, 1) = 0;
            end
        end
    end
    bw = MovingAverageThresh(tempim, floor(MATcrop/2), minArea, power);
    totalbw = bw;
    
    % Gets rid of the backgroudn noise.
    bw        = logical(bw);
    backavgIy = mean(double(Iy(~bw)));
    Iy        = uint16(floor(double(Iy) - backavgIy*ones(r, c)));
    backavgIc = mean(double(Ic(~bw)));
    Ic        = uint16(floor(double(Ic) - backavgIc*ones(r, c)));
    backavgIyy = mean(double(Iyy(~bw)));
    Iyy       = uint16(floor(double(Iyy) - backavgIyy*ones(r, c)));
    
    % This will record the total traces values of the cells.
    tottracks.yfp   = [tottracks.yfp mean2(double(Iy(bw)))];
    tottracks.cfp   = [tottracks.cfp mean2(double(Ic(bw)))];
    tottracks.ratio = [tottracks.ratio tottracks.cfp(i)/tottracks.yfp(i)];
    tottracks.rawy  = [tottracks.rawy mean2(double(Iyr(bw)))];
    tottracks.rawc  = [tottracks.rawc mean2(double(Icr(bw)))];
    tottracks.rawr  = [tottracks.rawr tottracks.rawc(i)/tottracks.rawy(i)];
    
    % We will now go through all the cells we have previously found.
    for j = 1:length(tracks)
        % We look at the cell only if it was found recently.
        maxt = max(tracks(j).frame);
        if maxt + maxf >= i
            % Find some important values.
            maxl = length(tracks(j).frame);
            a   = tracks(j).area(maxl);
            rad = sqrt(a/pi);
            cx  = tracks(j).x(maxl) + tracks(j).vx(maxl)*dt*(i - maxt);
            cy  = tracks(j).y(maxl) + tracks(j).vy(maxl)*dt*(i - maxt);
            if cx < 1
                cx = 1;
            elseif cx > c
                cx = c;
            end
            if cy <1
                cy = 1;
            elseif cy > r
                cy = r;
            end
            
            % Crop out the area around the cell in question.
            cropxl = int16(floor(cx - croprad*rad));
            if cropxl < 1
                cropxl = 1;
            end
            cropxh = int16(floor(cx + croprad*rad));
            if cropxh > c
                cropxh = c;
            end
            cropyl = int16(floor(cy - croprad*rad));
            if cropyl < 1
                cropyl = 1;
            end
            cropyh = int16(floor(cy + croprad*rad));
            if cropyh > r
                cropyh = r;
            end
            im_crop = im(cropyl:cropyh, cropxl:cropxh);
            im_crop = im_crop - min2(double(im_crop));
            im_crop = double(im_crop)/max2(double(im_crop));
            im_crop = uint16(floor((2^16-1)*im_crop));
            
            % Finding the binary image.
            bw       = MovingAverageThresh(im_crop, MATcrop, minArea, power);
            [r1, c1] = size(bw);
            tempim   = im_crop;
            for i1 = 1:r1
                for j1 = 1:c1
                    if bw(i1, j1) == 0
                        tempim(i1, j1) = 0;
                    end
                end
            end
            bw  = MovingAverageThresh(tempim, floor(MATcrop/2), minArea, power);
            bwt = zeros(r, c);
            bwt(cropyl:cropyh, cropxl:cropxh) = bw;
            bw  = bwt;
            
            % Finding the statistics of the cells in the cropped window.
            bw    = logical(bw);
            L     = bwlabel(bw, 4);
            stats = regionprops(L, 'Area', 'Centroid');
            
            % Finding the cell that shares coordinates of previous center.
            if ~isempty(stats)
                bwt = bwselect(bw, floor(cx), floor(cy), 4);
                if sum2(double(bwt)) ~= 0
                    bw = logical(bwt);
                    stats = regionprops(bw, 'Area', 'Centroid');
                else
                    counter = 0;
                    comp = (2*rad)^2;
                    for k = 1:length(stats)
                        comp1 = (stats(k).Centroid(1) - double(cx))^2 ...
                            + (stats(k).Centroid(2) - double(cy))^2;
                        if comp1 < comp
                            comp = comp1;
                            counter = k;
                        end
                    end
                    if counter ~= 0
                        bw = (L == counter);
                    else
                        bw = zeros(r, c);
                    end
                    stats = regionprops(bw, 'Area', 'Centroid');
                end
            end
            
            % Finding out if cell needs to be split.  Then splitting it.
            if ~isempty(stats)
                loop = (stats.Area - a > a*emarg) || ...
                    (pdist([stats.Centroid(1) stats.Centroid(2); cx cy]) > rad);
                if loop
                    loop = (pdist([stats.Centroid(1) stats.Centroid(2); cx cy]) > 1);
                end
                
                % If cell needs splitting, we crop the center of the cell.
                if loop
                    r1 = uint16(floor(stats.Centroid(2)));
                    c1 = uint16(floor(stats.Centroid(1)));
                    r1lower = r1 - bwcrppr;
                    if r1lower < 1
                        r1lower = 1;
                    end
                    r1upper = r1 + bwcrppr;
                    if r1upper > r
                        r1upper = r;
                    end
                    c1lower = c1 - bwcrppr;
                    if c1lower < 1
                        c1lower = 1;
                    end
                    c1upper = c1 + bwcrppr;
                    if c1upper > c
                        c1upper = c;
                    end
                    bw_crop = bw(r1lower:r1upper, c1lower:c1upper);
                    radius = 0;
                end
                bwt = bw;
                
                % We implement an infinite loop to seperate the cell.
                while loop
                    % Erosion to attempt to split the cells.
                    L = bwlabel(bw_crop, 4);
                    stats = regionprops(L, 'Centroid');
                    if ~isempty(stats)
                        radius = radius + 1;
                        circ = strel('disk', radius, 0);
                        bw_crop = imerode(bw_crop, circ);
                        
                        if ~sum2(double(bw_crop))
                            loop = 0;
                        end
                    else
                        loop = 0;
                    end
                    
                    % We will reassemble the cell and pick the center one.
                    bwt(r1lower:r1upper, c1lower:c1upper) = bw_crop;
                    bwt = bwareaopen(bwt, minArea, 4);
                    bwt = logical(bwt);
                    L = bwlabel(bwt, 4);
                    stats = regionprops(L, 'Centroid');
                    if ~isempty(stats)
                        counter = 0;
                        comp = Inf;
                        for k = 1:length(stats)
                            comp1 = (stats(k).Centroid(1) - double(cx))^2 + ...
                                (stats(k).Centroid(2) - double(cy))^2;
                            if comp1 < comp
                                comp = comp1;
                                counter = k;
                            end
                        end
                        
                        bwt = (L == counter);
                        
                        % We will now undo the erosion to check area.
                        bw_crop = bwt(r1lower:r1upper, c1lower:c1upper);
                        inv = 1 - bw_crop;
                        inv = imerode(inv, strel('disk', radius, 0));
                        bw_crop = 1 - inv;
                        
                       
                        bwt(r1lower:r1upper, c1lower:c1upper) = bw_crop;
                        
                        % We will now compute loop again.
                        L = bwlabel(bwt, 4);
                        stats = regionprops(L, 'Area', 'Centroid');
                        if loop
                            loop = (stats.Area - a > a*emarg) || ...
                                (pdist([stats.Centroid(1) stats.Centroid(2); cx cy]) > rad);
                        end
                        
                        if loop
                            bw_crop = bw(r1lower:r1upper, c1lower:c1upper);
                        end
                    else
                        loop = 0;
                    end
                end
                
                
                % Finding the stats of the single cell found again.
                bw = bwt;
                bw = logical(bw);
                L = bwlabel(bw, 4);
                stats = regionprops(L, 'Area', 'Centroid');
                stats12 = regionprops(imfill(L, 'holes'), 'Area');
            end
            
            
            
            % We can now record the data, IF it is near enough to previous.
            if ~isempty(stats)
                if (stats.Centroid(1) - double(cx))^2 + ...
                        (stats.Centroid(2) - double(cy))^2 < rad^2
                        bwt = logical(bw);
                        bw = imerode(bwt, strel('disk', erode, 0));
                    tracks(j).x     = [tracks(j).x stats.Centroid(1)];
                    tracks(j).y     = [tracks(j).y stats.Centroid(2)];
                    tracks(j).vx    = [tracks(j).vx (stats.Centroid(1) - tracks(j).x(maxl))/(dt*(i-maxt))];
                    tracks(j).vy    = [tracks(j).vy (stats.Centroid(2) - tracks(j).y(maxl))/(dt*(i-maxt))];
                    tracks(j).frame = [tracks(j).frame i];
                    tracks(j).area  = [tracks(j).area stats12.Area];
                    tracks(j).round = [tracks(j).round 4*pi*stats.Area/sum2(double(bwperim(bw,8)))];
                    tracks(j).yfp   = [tracks(j).yfp mean2(double(Iy(bw)))];
                    tracks(j).cfp   = [tracks(j).cfp mean2(double(Ic(bw)))];
                    tracks(j).yyfp  = [tracks(j).yyfp mean2(double(Iyy(bw)))];
                    tracks(j).ratio = [tracks(j).ratio tracks(j).cfp(maxl+1)/tracks(j).yfp(maxl+1)];
                    tracks(j).pratio = [tracks(j).pratio mean2(double(Ic(bw))./double(Iy(bw)))];
                    tracks(j).rawy  = [tracks(j).rawy mean2(double(Iyr(bw)))];
                    tracks(j).rawc  = [tracks(j).rawc mean2(double(Icr(bw)))];
                    tracks(j).rawyy = [tracks(j).rawyy mean2(double(Iyyr(bw)))];
                    tracks(j).rawr  = [tracks(j).rawr tracks(j).rawc(maxl+1)/tracks(j).rawy(maxl+1)];
                    tracks(j).prawr = [tracks(j).prawr mean2(double(Icr(bw))./double(Iyr(bw)))];
                    corrconstant    = 1; %((tracks(j).yyfp(1) - xtalk*tracks(j).cfp(1))/(tracks(j).yyfp(maxl+1) - xtalk*tracks(j).cfp(maxl+1)));
                    tracks(j).corr  = [tracks(j).corr corrconstant*(tracks(j).yfp(maxl+1) - alpha*tracks(j).yyfp(maxl+1) - beta*tracks(j).cfp(maxl+1))/(alpha*tracks(j).yyfp(maxl+1))];
                    tracks(j).pcorr = [tracks(j).pcorr corrconstant*mean2((double(Iy(bw)) - alpha.*double(Iyy(bw)) - beta.*double(Ic(bw)))./(alpha.*double(Iyy(bw))))];
                    tracks(j).rawcorr = [tracks(j).rawcorr corrconstant*(tracks(j).rawy(maxl+1) - alpha*tracks(j).rawyy(maxl+1) - beta*tracks(j).rawc(maxl+1))/(alpha*tracks(j).rawyy(maxl+1))];
                    tracks(j).prawcorr = [tracks(j).prawcorr corrconstant*mean2((double(Iyr(bw)) - alpha.*double(Iyyr(bw)) - beta.*double(Icr(bw)))./(alpha.*double(Iyyr(bw))))];
                    masks{j} = bw;
                    totalbw(bwt) = 0;
                end
            end
        end
    end
    
    % Now we will search for cells that we have not yet found.  New cells.
    totalbw = bwareaopen(totalbw, minArea, 4);
    L = bwlabel(totalbw, 4);
    stats = regionprops(L, 'Area', 'Centroid');
    if ~isempty(stats)
        % Runs through all the found cells to get rid of doubles and
        % connect if necessary.
        for k = 1:length(tracks)
            maxt = max(tracks(k).frame);
            % This will get rid of the cells that are dubles.
            if maxt == i
                maxl = length(tracks(k).frame);
                a = tracks(k).area(maxl);
                rad = sqrt(a/pi);
                for j = 1:length(stats)
                    dist = pdist([stats(j).Centroid(1) stats(j).Centroid(2); ...
                        tracks(k).x(maxl) tracks(k).y(maxl)]);
                    if dist < rad*2
                        totalbw(L == j) = 0;
                        stats(j).Area = Inf;
                    end
                end
            % This will cennect the cells that are the same.
            elseif maxt + maxf2 >= i
                maxl = length(tracks(k).frame);
                a = tracks(k).area(maxl);
                rad = sqrt(a/pi);
                cx  = tracks(k).x(maxl) + tracks(k).vx(maxl)*dt*(i - maxt);
                cy  = tracks(k).y(maxl) + tracks(k).vy(maxl)*dt*(i - maxt);
                if cx < 1
                    cx = 1;
                elseif cx > c
                    cx = c;
                end
                if cy <1
                    cy = 1;
                elseif cy > r
                    cy = r;
                end
                
                for j = 1:length(stats)
                    if (stats(j).Area - a < a*emarg) && ...
                            ((stats(j).Centroid(1) - cx)^2 + (stats(j).Centroid(2) - cy)^2 < rad^2)
                        bw = (L == j);
                        bw = logical(bw);
                        stats12 = regionprops(imfill(bw, 'holes'), 'Area');
                        bwt = bw;
                        bw = imerode(bwt, strel('disk', erode, 0));
                        tracks(k).x     = [tracks(k).x stats(j).Centroid(1)];
                        tracks(k).y     = [tracks(k).y stats(j).Centroid(2)];
                        tracks(k).vx    = [tracks(k).vx (stats(j).Centroid(1) - tracks(k).x(maxl))/(dt*(i-maxt))];
                        tracks(k).vy    = [tracks(k).vy (stats(j).Centroid(2) - tracks(k).y(maxl))/(dt*(i-maxt))];
                        tracks(k).frame = [tracks(k).frame i];
                        tracks(k).area  = [tracks(k).area stats12.Area];
                        tracks(k).round = [tracks(k).round 4*pi*stats12.Area/sum2(double(bwperim(bw, 8)))];
                        tracks(k).yfp   = [tracks(k).yfp mean2(double(Iy(bw)))];
                        tracks(k).cfp   = [tracks(k).cfp mean2(double(Ic(bw)))];
                        tracks(k).yyfp  = [tracks(k).yyfp mean2(double(Iyy(bw)))];
                        tracks(k).ratio = [tracks(k).ratio tracks(k).cfp(maxl+1)/tracks(k).yfp(maxl+1)];
                        tracks(k).pratio = [tracks(k).pratio mean2(double(Ic(bw))./double(Iy(bw)))];
                        tracks(k).rawy  = [tracks(k).rawy mean2(double(Iyr(bw)))];
                        tracks(k).rawc  = [tracks(k).rawc mean2(double(Icr(bw)))];
                        tracks(k).rawyy = [tracks(k).rawyy mean2(double(Iyyr(bw)))];
                        tracks(k).rawr  = [tracks(k).rawr tracks(k).rawc(maxl+1)/tracks(k).rawy(maxl+1)];
                        tracks(k).prawr = [tracks(k).prawr mean2(double(Icr(bw))./double(Iyr(bw)))];
                        corrconstant    = 1; %((tracks(k).yyfp(1) - xtalk*tracks(k).cfp(1))/(tracks(k).yyfp(maxl+1) - xtalk*tracks(k).cfp(maxl+1)));
                        tracks(k).corr  = [tracks(k).corr corrconstant*(tracks(k).yfp(maxl+1) - alpha*tracks(k).yyfp(maxl+1) - beta*tracks(k).cfp(maxl+1))/(alpha*tracks(k).yyfp(maxl+1))];
                        tracks(k).pcorr = [tracks(k).pcorr corrconstant*mean2((double(Iy(bw)) - alpha.*double(Iyy(bw)) - beta.*double(Ic(bw)))./(alpha.*double(Iyy(bw))))];
                        tracks(k).rawcorr = [tracks(k).rawcorr corrconstant*(tracks(k).rawy(maxl+1) - alpha*tracks(k).rawyy(maxl+1) - beta*tracks(k).rawc(maxl+1))/(alpha*tracks(k).rawyy(maxl+1))];
                        tracks(k).prawcorr = [tracks(k).prawcorr corrconstant*mean2((double(Iyr(bw)) - alpha.*double(Iyyr(bw)) - beta.*double(Icr(bw)))./(alpha.*double(Iyyr(bw))))];
                        masks{k} = bw;
                        total(bwt) = 0;
                    end
                end
            end
        end
    end
    
    % Now we enter into the section where we add the new cells.
    L = bwlabel(totalbw);
    stats = regionprops(L, 'Area', 'Centroid');
    temptracks = length(tracks);
    if ~isempty(stats)
        for j = 1:length(stats)
            bw      = (L == j);
            bw      = logical(bw);
            stats12 = regionprops(imfill(bw, 'holes'), 'Area');
            bwt = bw;
            bw = imerode(bwt, strel('disk', erode, 0));
            tracks(temptracks + j).x     = stats(j).Centroid(1);
            tracks(temptracks + j).y     = stats(j).Centroid(2);
            tracks(temptracks + j).vx    = 0;
            tracks(temptracks + j).vy    = 0;
            tracks(temptracks + j).frame = i;
            tracks(temptracks + j).area  = stats12.Area;
            tracks(temptracks + j).round = 4*pi*stats12.Area/sum2(double(bwperim(bw, 8)));
            tracks(temptracks + j).yfp   = mean2(double(Iy(bw)));
            tracks(temptracks + j).cfp   = mean2(double(Ic(bw)));
            tracks(temptracks + j).yyfp  = mean2(double(Iyy(bw)));
            tracks(temptracks + j).ratio = tracks(temptracks+j).cfp/tracks(temptracks+j).yfp;
            tracks(temptracks + j).pratio = mean2(double(Ic(bw))./double(Iy(bw)));
            tracks(temptracks + j).rawy  = mean2(double(Iyr(bw)));
            tracks(temptracks + j).rawc  = mean2(double(Icr(bw)));
            tracks(temptracks + j).rawyy = mean2(double(Iyyr(bw)));
            tracks(temptracks + j).rawr  = tracks(temptracks+j).rawc/tracks(temptracks+j).rawy;
            tracks(temptracks + j).prawr = mean2(double(Icr(bw))./double(Iyr(bw)));
            tracks(temptracks + j).corr = (tracks(temptracks + j).yfp - alpha*tracks(temptracks + j).yyfp - beta*tracks(temptracks + j).cfp)/(alpha*tracks(temptracks + j).yyfp);
            tracks(temptracks + j).pcorr = mean((double(Iy(bw)) - alpha.*double(Iyy(bw)) - beta.*double(Ic(bw)))/(alpha.*double(Iyy(bw))));
            tracks(temptracks + j).rawcorr = (tracks(temptracks + j).rawy - alpha*tracks(temptracks + j).rawyy - beta*tracks(temptracks + j).rawc)/(alpha*tracks(temptracks + j).rawyy);
            tracks(temptracks + j).prawcorr = mean((double(Iyr(bw)) - alpha.*double(Iyyr(bw)) - beta.*double(Icr(bw)))/(alpha.*double(Iyyr(bw))));
            masks{temptracks + j}        = bw;
        end
    end
    
    % We will now add a frame to the movie.
    im = double(im)/max2(double(im));
    im = repmat(im, [1 1 3]);
    totalmasks = zeros(r, c);
    totalmasks = repmat(totalmasks, [1 1 3]);
    for j = 1:length(tracks)
        if sum(tracks(j).frame == i)
            for t = 1:3
                if mod(j, modder) ~= 0
                    totalmasks(:, :, t) = totalmasks(:, :, t) + cmap(mod(j, modder), t)*masks{j};
                else
                    totalmasks(:, :, t) = totalmasks(:, :, t) + cmap(modder, t)*masks{j};
                end
            end
        end
    end
    aviobj = addframe(aviobj, [im totalmasks]);
end
aviobj = close(aviobj);
save([imdir, 'tracks.mat'], 'tracks');
clearvars -except rootdir yy_name y_name c_name dt startImg numImg tau tau2 positions alpha beta; clf; clc;
 end
%--------------------------------------------------------------------------


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
% Modified by Allyson Sgro                                       2012-09-17
% Modified by Allyson Sgro                                       2012-10-08
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

