% AUTOCONTRAST  Automatically adjusts contrast of images to optimum level.
%    e.g. autocontrast('Sunset.jpg','Output.jpg')
% maxValue is 255 for 8-bit and 65535 for 16 bit image
function output_img = autocontrast(img,numBits)

low_limit=0.0001;
up_limit=0.9999;
[m1 n1 r1]=size(img);
img=double(img);
%--------------------calculation of vmin and vmax----------------------
for k=1:r1
    arr=sort(reshape(img(:,:,k),m1*n1,1));
    v_min(k)=arr(ceil(low_limit*m1*n1));
    v_max(k)=arr(ceil(up_limit*m1*n1));
end
%----------------------------------------------------------------------
if r1==3
    v_min=rgb2ntsc(v_min);
    v_max=rgb2ntsc(v_max);
end
%----------------------------------------------------------------------
tmpimg=(img-v_min(1))/(v_max(1)-v_min(1));
if numBits==8
    output_img = uint8(tmpimg.*255);

elseif numBits==16
    output_img = uint16(tmpimg.*65535);
else
    output_img = img;    
end
