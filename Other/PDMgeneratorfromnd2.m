clear all; clc;
tic
cd(userpath);

inDir = 'C:\Users\dplad\Desktop\testfolder\*.nd2';
Folder = 'C:\Users\dplad\Desktop\testfolder\';

srcFiles = dir(inDir);
PercentileThresh(1,1:4) = [99 99 99 99];
First = 1; %%This and the next variable indicate which channels to use for PDM.
Second = 2;
Channels = 2; %%The total number of channels in the image.
bgsubstrel = strel('disk',50);
SAVE = 1; %Will save the final PDM images  to folder if SAVE=1. If SAVE=0, no images will be saved.

Ch1_translation = [0,0,0];
Ch2_translation = [0,0,0];
Ch3_translation = [0,0,0];
Ch4_translation = [0,0,0];



for f = 2:length(srcFiles)
    progress = (f/length(srcFiles))*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    disp(progress2)
    
    filename = strcat(Folder,srcFiles(f).name);
    
    cd(Folder);
    
    newF = strcat(srcFiles(f).name(1:(end-4)));
    if SAVE>0, mkdir(newF); else end;
    loop = num2str(f);
    slash = '\';
    if SAVE>0, ImageDir = strcat(Folder,newF,slash); else end;
    
I = bfopen(filename);

Res = length(I{1,1}{1,1});
Slices = (length(I{1,1})/Channels);
Blank = zeros(Res,Res,Slices);
Planes = Channels*Slices;

Ch1 = uint16(Blank);
Ch2 = uint16(Blank);
Ch3 = uint16(Blank);
Ch4 = uint16(Blank);
ZeroMask = double(Blank);

First_DM = uint16(Blank);
Second_DM = uint16(Blank);

PDM_pos = uint16(Blank);
PDM_neg = uint16(Blank);

for i = 1:Slices
    Ch1_planes(i,1) = 1+(Channels*i-Channels);
    Ch2_planes(i,1) = 2+(Channels*i-Channels);
    if Channels==3, Ch3_planes(i,1) = 3+(Channels*i-Channels); else end;
    if Channels==4, Ch4_planes(i,1) = 4+(Channels*i-Channels); else end;
end

for m = 1:Slices
    Ch1(:,:,m) = I{1,1}{Ch1_planes(m,1),1};
    Ch2(:,:,m) = I{1,1}{Ch2_planes(m,1),1};
    if Channels==3, Ch3(:,:,m) = I{1,1}{Ch3_planes(m,1),1}; else end
    if Channels==4, Ch4(:,:,m) = I{1,1}{Ch4_planes(m,1),1}; else end
end

% Ch1 = rescale(Ch1); Ch1_bg = imopen(Ch1,bgsubstrel); Ch1 = Ch1 - Ch1_bg;
% Ch2 = rescale(Ch2); Ch2_bg = imopen(Ch2,bgsubstrel); Ch2 = Ch2 - Ch2_bg;
% if Channels==3, Ch3_bg = imopen(Ch3,bgsubstrel); Ch3 = im2double(Ch3 - Ch3_bg); else end
% if Channels==4, Ch4_bg = imopen(Ch4,bgsubstrel); Ch4 = im2double(Ch4 - Ch4_bg); else end

Ch1_thresh = prctile(Ch1,PercentileThresh(1,1),'all');
maskch1 = Ch1>Ch1_thresh; Ch1_seg = double(Ch1).*maskch1;
Ch1_seg = imtranslate(Ch1_seg,Ch1_translation);

Ch2_thresh = prctile(Ch2,PercentileThresh(1,2),'all');
maskch2 = Ch2>Ch2_thresh; Ch2_seg = double(Ch2).*maskch2;
Ch2_seg = imtranslate(Ch2_seg,Ch2_translation);

if Channels==3,
Ch3_thresh = prctile(Ch3,PercentileThresh(1,3),'all');
maskch3 = Ch3>Ch3_thresh; Ch3_seg = Ch3.*maskch3;
Ch3_seg = imtranslate(Ch3_seg,Ch3_translation);
else end

if Channels==4,
Ch4_thresh = prctile(Ch4,PercentileThresh(1,4),'all');
maskch4 = Ch4>Ch4_thresh; Ch4_seg = Ch4.*maskch4;
Ch4_seg = imtranslate(Ch4_seg,Ch4_translation);
else end

if First==1, FirstCh = Ch1_seg;
elseif First==2, FirstCh = Ch2_seg;
elseif First==3, FirstCh = Ch3_seg;
elseif First==4, FirstCh = Ch4_seg;
else end

if Second==1, SecondCh = Ch1_seg;
elseif Second==2, SecondCh = Ch2_seg;
elseif Second==3, SecondCh = Ch3_seg;
elseif Second==4, SecondCh = Ch4_seg;
else end

FirstCh_Mean = mean(FirstCh(:)>0,'all');
SecondCh_Mean = mean(SecondCh(:)>0,'all');

for y = 1:Res
    for x = 1:Res
        for z = 1:Slices
        First_DM(y,x,z) = FirstCh(y,x,z) - FirstCh_Mean;
        Second_DM(y,x,z) = SecondCh(y,x,z) - SecondCh_Mean;
        PDM(y,x,z) = First_DM(y,x,z)*Second_DM(y,x,z);
        
        if PDM(y,x,z) >0, PDM_pos(y,x,z) = PDM(y,x,z);
        else PDM_pos(y,x,z) == 0;
        end
        if PDM(y,x,z) <0, PDM_neg(y,x,z) = PDM(y,x,z);
        else PDM_neg(y,x,z) == 0;
        end
        end
    end
end

for y = 1:Res
    for x = 1:Res
        for z = 1:Slices
            if First_DM(y,x,z) <0 || Second_DM(y,x,z) <0, ZeroMask(y,x,z) = 0;
            else ZeroMask(y,x,z) = 1;
            end
        end
    end
end

ZeroMask = logical(ZeroMask);
PDM_pos2 = double(PDM_pos).*ZeroMask;
PDM_neg2 = double(-PDM_neg).*ZeroMask;

%PDM_pos16 = im2uint16(PDM_pos2);
%PDM_neg16 = im2uint16(-PDM_neg2);
%PDM_neg = -PDM_neg;

if SAVE>0,
    cd(ImageDir);
    bfsave(PDM_pos2,'PDM_pos.ome.tiff','dimensionOrder','XYZTC');
    bfsave(PDM_neg2,'PDM_neg.ome.tiff','dimensionOrder','XYZTC');
else end
end

toc