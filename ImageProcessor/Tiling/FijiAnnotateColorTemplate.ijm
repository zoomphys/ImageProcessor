//////////////////////////////////
// Author: Zoom Nguyen
// Email: duongnh@gmail.com
// Updated: 2012/12/01, 2013/05/03
//////////////////////////////////

// set variables from external parameters
scale_known_um = $umScale;
framePerMin = $frmPerMin;

//////////////////////////////////
scalebar_width = 100;
time_locx = 20;
time_locy = 50;

scalebar_loc = "Lower Right";
scale_distance_pixel = 1;

// graylevel_time = 0;
// color_scalebar = "Black";
graylevel_time = 255;
color_scalebar = "White";

// row in results table containing slice number
slice_row = 0;

slice = getResult("var", slice_row);
time = floor((slice-1)/framePerMin);

run("Duplicate...","title=Duplicated");

setFont("SansSerif", 24, "bold antialiased");
setForegroundColor(graylevel_time, graylevel_time, graylevel_time);

///////////////////////////
run("Gaussian Blur...", "sigma=1");
run("Fire");
run("RGB Color");

drawString(time+" min", time_locx, time_locy);
run("Colors...", "foreground=white background=black selection=yellow");

run("Set Scale...", "distance="+scale_distance_pixel+" known="+scale_known_um+" pixel=1 unit=um");
run("Scale Bar...", "width="+scalebar_width+" height=8 font=24 color="+color_scalebar+" background=None location=["+scalebar_loc+"] bold");

// Don't select anything so that the whole image can be saved in the next step
run("Select None");