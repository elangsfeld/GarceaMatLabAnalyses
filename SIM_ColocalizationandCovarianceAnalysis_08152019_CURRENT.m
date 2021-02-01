 clear all; clc;
tic
cd(userpath);

inDir = 'C:\Users\dope8688\Desktop\fish edu\*.nd2';
Folder = 'C:\Users\dope8688\Desktop\fish edu\';

srcFiles = dir(inDir);%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LPercentileThresh(1,1:4) = [99.75 99.75 99.75 99.75];
UPercentileThresh(1,1:4) = [100 100 100 100];
LStaticThresh(1,1:4) = [5000 4000 5000 1000];
UStaticThresh(1,1:4) = [60000 60000 60000 60000];
Channels = 4;
%bgsubstrel = strel('disk',50);
SAVE_IMAGE = 0; %Will save separated and thresholded images to folder if =1. If =0, no images will be saved (speeds up analysis).
SAVE_PCC = 0; % Will save 2-column matrices for each PCC comparison (1=yes, 0=no).
THRESH_TYPE = 0; %Determines which thresholding method is used. 0=Percentile, 1=Static Value.
%%%%XYZ TRANSLATION VECTOR%%%%%%%%
Ch1_translation = [0,0,0];
Ch2_translation = [0,0,0];
Ch3_translation = [0,0,0];
Ch4_translation = [0,0,0];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
colNames = {'File','PCC12','PCC13','PCC14','PCC23','PCC24','PCC34','M1o12','M1o13','M1o14','M1o23','M1o24','M1o34','M2o12','M2o13','M2o14','M2o23','M2o24','M2o34','PCCPIX1o2','PCCPIX1o3','PCCPIX1o4','PCCPIX2o3','PCCPIX2o4','PCCPIX3o4','MCCPIX1o2','MCCPIX1o3','MCCPIX1o4','MCCPIX2o3','MCCPIX2o4','MCCPIX3o4','PIXCh1','PIXCh2','PIXCh3','PIXCh4'};
colNames2 = colNames';
varTypes = {'string','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};
PCCMOC = table('Size',[length(srcFiles) 35],'VariableNames',colNames,'VariableTypes',varTypes);
%%%

for f = 1:length(srcFiles)
progress = (f/length(srcFiles))*100;
progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
disp(progress2)
    
    filename = strcat(Folder,srcFiles(f).name);

    cd(Folder);
    
    newF = strcat(srcFiles(f).name(1:(end-4)));
    PCCMOC.File(f,1) = newF;
    if SAVE_IMAGE > 0 || SAVE_PCC > 0, mkdir(newF); else end;
    loop = num2str(f);
    slash = '\';
    if SAVE_IMAGE > 0 || SAVE_PCC > 0, ImageDir = strcat(Folder,newF,slash); else end;
    
I = bfopen(filename);

Res = length(I{1,1}{1,1});
Slices = (length(I{1,1})/Channels);  %%Calculates number of Z positions
Blank = zeros(Res,Res,Slices);
Planes = Channels*Slices; %%Back-calculates number of elements in "I".

%%Blanks for later%%
Ch1 = uint16(Blank);
Ch2 = uint16(Blank);
Ch3 = uint16(Blank);
Ch4 = uint16(Blank);
%%

%%%%Creates arrays that index positions in "I" that correspond to each fluorescence channel.
for i = 1:Slices
    Ch1_planes(i,1) = 1+(Channels*i-Channels);
    Ch2_planes(i,1) = 2+(Channels*i-Channels);
    if Channels>2, Ch3_planes(i,1) = 3+(Channels*i-Channels); else end
    if Channels>3, Ch4_planes(i,1) = 4+(Channels*i-Channels); else end
end

%%%%Separates each fluorescence channel into its own stack of images.
for m = 1:Slices
    Ch1(:,:,m) = I{1,1}{Ch1_planes(m,1),1};
    Ch2(:,:,m) = I{1,1}{Ch2_planes(m,1),1};
    if Channels>2, Ch3(:,:,m) = I{1,1}{Ch3_planes(m,1),1}; else end
    if Channels>3, Ch4(:,:,m) = I{1,1}{Ch4_planes(m,1),1}; else end
end

%Ch1_bg = imopen(Ch1,bgsubstrel); Ch1 = Ch1 - Ch1_bg;
%Ch2_bg = imopen(Ch2,bgsubstrel); Ch2 = Ch2 - Ch2_bg;
%if Channels>2, Ch3_bg = imopen(Ch3,bgsubstrel); Ch3 = Ch3 - Ch3_bg; else end
%if Channels>3, Ch4_bg = imopen(Ch4,bgsubstrel); Ch4 = Ch4 - Ch4_bg; else end


%%XYTZC%%

if SAVE_IMAGE > 0, cd(ImageDir); else end;

%%%%These save the raw stack image of each channel. These statements also 
%%%%create a mask image for each channel based on the percentile threshold 
%%%%indicated at the top, then saves that masked/segmented image.
    if SAVE_IMAGE > 0, bfsave(Ch1,'Ch1-stack.ome.tiff','dimensionOrder','XYZTC'); else end;
    
    if THRESH_TYPE == 0, 
        LCh1_thresh = prctile(Ch1,LPercentileThresh(1,1),'all'); 
        UCh1_thresh = prctile(Ch1,UPercentileThresh(1,1),'all');
    else
        LCh1_thresh = LStaticThresh(1,1);
        UCh1_thresh = UStaticThresh(1,1);
    end
    maskch1 = (UCh1_thresh>Ch1) & (Ch1>LCh1_thresh); Ch1_seg = uint16(double(Ch1).*maskch1);
    Ch1_seg = imtranslate(Ch1_seg,Ch1_translation);
    if SAVE_IMAGE > 0, bfsave(Ch1_seg,'Ch1segmented-stack.ome.tiff','dimensionOrder','XYZTC'); else end;

    if SAVE_IMAGE > 0, bfsave(Ch2,'Ch2-stack.ome.tiff','dimensionOrder','XYZTC'); else end;
    
    if THRESH_TYPE == 0, 
        LCh2_thresh = prctile(Ch2,LPercentileThresh(1,2),'all');
        UCh2_thresh = prctile(Ch2,UPercentileThresh(1,2),'all');
    else
        LCh2_thresh = LStaticThresh(1,2);
        UCh2_thresh = UStaticThresh(1,2);
    end
    maskch2 = (UCh2_thresh>Ch2) & (Ch2>LCh2_thresh); Ch2_seg = uint16(double(Ch2).*maskch2);
    Ch2_seg = imtranslate(Ch2_seg,Ch2_translation);
    if SAVE_IMAGE > 0, bfsave(Ch2_seg,'Ch2segmented-stack.ome.tiff','dimensionOrder','XYZTC'); else end;

if Channels>2,
    if SAVE_IMAGE > 0, bfsave(Ch3,'Ch3-stack.ome.tiff','dimensionOrder','XYZTC'); else end;
    
    if THRESH_TYPE == 0, 
        LCh3_thresh = prctile(Ch3,LPercentileThresh(1,3),'all');
        UCh3_thresh = prctile(Ch3,UPercentileThresh(1,3),'all');
    else
        LCh3_thresh = LStaticThresh(1,3);
        UCh3_thresh = UStaticThresh(1,4);
    end
    maskch3 = (UCh3_thresh>Ch3) & (Ch3>LCh3_thresh); Ch3_seg = uint16(double(Ch3).*maskch3);
    Ch3_seg = imtranslate(Ch3_seg,Ch3_translation);
    if SAVE_IMAGE > 0, bfsave(Ch3_seg,'Ch3segmented-stack.ome.tiff','dimensionOrder','XYZTC'); else end;
else end

if Channels>3,
    if SAVE_IMAGE > 0, bfsave(Ch4,'Ch4-stack.ome.tiff','dimensionOrder','XYZTC'); else end;
    
    if THRESH_TYPE == 0, 
        LCh4_thresh = prctile(Ch4,LPercentileThresh(1,4),'all');
        UCh4_thresh = prctile(Ch4,UPercentileThresh(1,4),'all');
    else
        LCh4_thresh = LStaticThresh(1,4);
        UCh4_thresh = UStaticThresh(1,4);
    end
    maskch4 = (UCh4_thresh>Ch4) & (Ch4>LCh4_thresh); Ch4_seg = uint16(double(Ch4).*maskch4);
    Ch4_seg = imtranslate(Ch4_seg,Ch4_translation);
    if SAVE_IMAGE > 0, bfsave(Ch4_seg,'Ch4segmented-stack.ome.tiff','dimensionOrder','XYZTC'); else end;
else end

%%%%These linearize each of the segmented images into a single column vector.
ReshapeL = numel(Ch1_seg); ReshapeM = zeros(ReshapeL,2);

Ch1_re = reshape(Ch1_seg,[ReshapeL 1]);
Ch2_re = reshape(Ch2_seg,[ReshapeL 1]);
if Channels>2, Ch3_re = reshape(Ch3_seg,[ReshapeL 1]); else end
if Channels>3, Ch4_re = reshape(Ch4_seg,[ReshapeL 1]); else end

%%%%These calculate number of above-threshold pixels in each channel.
if sum(Ch1_re(:))>0, PCCMOC(f,32) = table(nnz(Ch1_re(:))); else end
if sum(Ch2_re(:))>0, PCCMOC(f,33) = table(nnz(Ch2_re(:))); else end
if Channels>2, if sum(Ch3_re(:))>0, PCCMOC(f,34) = table(nnz(Ch3_re(:))); else end
else end
if Channels>3, if sum(Ch4_re(:))>0, PCCMOC(f,35) = table(nnz(Ch4_re(:))); else end
else end

%%%%These pre-allocate arrays for the upcoming 0,0 pair screening.
Ch1o2_P = ReshapeM;
Ch1o2_M = ReshapeM;
if Channels>2,
    Ch1o3_P = ReshapeM; Ch2o3_P = ReshapeM;
    Ch1o3_M = ReshapeM; Ch2o3_M = ReshapeM;
else end
if Channels>3,
    Ch1o4_P = ReshapeM; Ch2o4_P = ReshapeM; Ch3o4_P = ReshapeM;
    Ch1o4_M = ReshapeM; Ch2o4_M = ReshapeM; Ch3o4_M = ReshapeM;
else end

%%%%This loop identifies rows where EITHER channel has 0 value (for PCC).
parfor w = 1:ReshapeL
    if (Ch1_re(w,1) == 0) || (Ch2_re(w,1) == 0), Ch1o2_P(w,:) = [NaN NaN];
    else Ch1o2_P(w,:) = [Ch1_re(w,1) Ch2_re(w,1)];
    end
    if Channels>2,
        if (Ch1_re(w,1) == 0) || (Ch3_re(w,1) == 0), Ch1o3_P(w,:) = [NaN NaN];
        else Ch1o3_P(w,:) = [Ch1_re(w,1) Ch3_re(w,1)];
        end
        if (Ch2_re(w,1) == 0) || (Ch3_re(w,1) == 0), Ch2o3_P(w,:) = [NaN NaN];
        else Ch2o3_P(w,:) = [Ch2_re(w,1) Ch3_re(w,1)]; 
        end
    else
    end
    if Channels>3
        if (Ch1_re(w,1) == 0) || (Ch4_re(w,1) == 0), Ch1o4_P(w,:) = [NaN NaN];
        else Ch1o4_P(w,:) = [Ch1_re(w,1) Ch4_re(w,1)];
        end
        if (Ch2_re(w,1) == 0) || (Ch4_re(w,1) == 0), Ch2o4_P(w,:) = [NaN NaN];
        else Ch2o4_P(w,:) = [Ch2_re(w,1) Ch4_re(w,1)];
        end
        if (Ch3_re(w,1) == 0) || (Ch4_re(w,1) == 0), Ch3o4_P(w,:) = [NaN NaN];
        else Ch3o4_P(w,:) = [Ch3_re(w,1) Ch4_re(w,1)];
        end
    else
    end
end

%%%%This loop identifies rows where BOTH channels have 0 value (for MOC)
parfor w = 1:ReshapeL
    if (Ch1_re(w,1) == 0) && (Ch2_re(w,1) == 0), Ch1o2_M(w,:) = [NaN NaN];
    else Ch1o2_M(w,:) = [Ch1_re(w,1) Ch2_re(w,1)];
    end
    if Channels>2,
        if (Ch1_re(w,1) == 0) && (Ch3_re(w,1) == 0), Ch1o3_M(w,:) = [NaN NaN];
        else Ch1o3_M(w,:) = [Ch1_re(w,1) Ch3_re(w,1)];
        end
        if (Ch2_re(w,1) == 0) && (Ch3_re(w,1) == 0), Ch2o3_M(w,:) = [NaN NaN];
        else Ch2o3_M(w,:) = [Ch2_re(w,1) Ch3_re(w,1)];
        end
    else
    end
    if Channels>3,
        if (Ch1_re(w,1) == 0) && (Ch4_re(w,1) == 0), Ch1o4_M(w,:) = [NaN NaN];
        else Ch1o4_M(w,:) = [Ch1_re(w,1) Ch4_re(w,1)];
        end
        if (Ch2_re(w,1) == 0) && (Ch4_re(w,1) == 0), Ch2o4_M(w,:) = [NaN NaN];
        else Ch2o4_M(w,:) = [Ch2_re(w,1) Ch4_re(w,1)];
        end
        if (Ch3_re(w,1) == 0) && (Ch4_re(w,1) == 0), Ch3o4_M(w,:) = [NaN NaN];
        else Ch3o4_M(w,:) = [Ch3_re(w,1) Ch4_re(w,1)];
        end
    else
    end
end

%%%%These calculate the number of pixels being used for PCC and MCC analyses.
PCCMOC(f,20) = table(ReshapeL - sum(isnan(Ch1o2_P(:,1))));
PCCMOC(f,26) = table(ReshapeL - sum(isnan(Ch1o2_M(:,1))));
if Channels>2,
PCCMOC(f,21) = table(ReshapeL - sum(isnan(Ch1o3_P(:,1))));
PCCMOC(f,23) = table(ReshapeL - sum(isnan(Ch2o3_P(:,1))));
PCCMOC(f,27) = table(ReshapeL - sum(isnan(Ch1o3_M(:,1))));
PCCMOC(f,29) = table(ReshapeL - sum(isnan(Ch2o3_M(:,1))));
else end
if Channels>3,
    PCCMOC(f,22) = table(ReshapeL - sum(isnan(Ch1o4_P(:,1))));
    PCCMOC(f,24) = table(ReshapeL - sum(isnan(Ch2o4_P(:,1))));
    PCCMOC(f,25) = table(ReshapeL - sum(isnan(Ch3o4_P(:,1))));
    PCCMOC(f,28) = table(ReshapeL - sum(isnan(Ch1o4_M(:,1))));
    PCCMOC(f,30) = table(ReshapeL - sum(isnan(Ch2o4_M(:,1))));
    PCCMOC(f,31) = table(ReshapeL - sum(isnan(Ch3o4_M(:,1))));
else end

%%%%This remove NaN values in each comparison for PCC and MCC analyses.
Ch1o2_P = rmmissing(Ch1o2_P); Ch1o2_M = rmmissing(Ch1o2_M);

if Channels>2,
    Ch1o3_P = rmmissing(Ch1o3_P); Ch1o3_M = rmmissing(Ch1o3_M);
    Ch2o3_P = rmmissing(Ch2o3_P); Ch2o3_M = rmmissing(Ch2o3_M);
else end
if Channels>3,
    Ch1o4_P = rmmissing(Ch1o4_P); Ch1o4_M = rmmissing(Ch1o4_M);
    Ch2o4_P = rmmissing(Ch2o4_P); Ch2o4_M = rmmissing(Ch2o4_M);
    Ch3o4_P = rmmissing(Ch3o4_P); Ch3o4_M = rmmissing(Ch3o4_M);
else end

%%%%This calculates M1 and M2 Mander's Overlap Coefficients.
Ch1o2_mask = zeros(length(Ch1o2_M),2);
Ch1o2_mask(:,1) = Ch1o2_M(:,1)>0; Ch1o2_mask(:,2) = Ch1o2_M(:,2)>0;
PCCMOC(f,8) = table(sum((Ch1o2_M(:,1).*Ch1o2_mask(:,2))/sum(Ch1o2_M(:,1))));
PCCMOC(f,14) = table(sum((Ch1o2_M(:,2).*Ch1o2_mask(:,1))/sum(Ch1o2_M(:,2))));

if Channels>2,
    Ch1o3_mask = zeros(length(Ch1o3_M),2);
    Ch1o3_mask(:,1) = Ch1o3_M(:,1)>0; Ch1o3_mask(:,2) = Ch1o3_M(:,2)>0;
    PCCMOC(f,9) = table(sum((Ch1o3_M(:,1).*Ch1o3_mask(:,2))/sum(Ch1o3_M(:,1))));
    PCCMOC(f,15) = table(sum((Ch1o3_M(:,2).*Ch1o3_mask(:,1))/sum(Ch1o3_M(:,2))));
    Ch2o3_mask = zeros(length(Ch2o3_M),2);
    Ch2o3_mask(:,1) = Ch2o3_M(:,1)>0; Ch2o3_mask(:,2) = Ch2o3_M(:,2)>0;
    PCCMOC(f,11) = table(sum((Ch2o3_M(:,1).*Ch2o3_mask(:,2))/sum(Ch2o3_M(:,1))));
    PCCMOC(f,17) = table(sum((Ch2o3_M(:,2).*Ch2o3_mask(:,1))/sum(Ch2o3_M(:,2))));
else end
    
if Channels>3,
    Ch1o4_mask = zeros(length(Ch1o4_M),2);
    Ch1o4_mask(:,1) = Ch1o4_M(:,1)>0; Ch1o4_mask(:,2) = Ch1o4_M(:,2)>0;
    PCCMOC(f,10) = table(sum((Ch1o4_M(:,1).*Ch1o4_mask(:,2))/sum(Ch1o4_M(:,1))));
    PCCMOC(f,16) = table(sum((Ch1o4_M(:,2).*Ch1o4_mask(:,1))/sum(Ch1o4_M(:,2))));
    Ch2o4_mask = zeros(length(Ch2o4_M),2);
    Ch2o4_mask(:,1) = Ch2o4_M(:,1)>0; Ch2o4_mask(:,2) = Ch2o4_M(:,2)>0;
    PCCMOC(f,12) = table(sum((Ch2o4_M(:,1).*Ch2o4_mask(:,2))/sum(Ch2o4_M(:,1))));
    PCCMOC(f,18) = table(sum((Ch2o4_M(:,2).*Ch2o4_mask(:,1))/sum(Ch2o4_M(:,2))));
    Ch3o4_mask = zeros(length(Ch3o4_M),2);
    Ch3o4_mask(:,1) = Ch3o4_M(:,1)>0; Ch3o4_mask(:,2) = Ch3o4_M(:,2)>0;
    PCCMOC(f,13) = table(sum((Ch3o4_M(:,1).*Ch3o4_mask(:,2))/sum(Ch3o4_M(:,1))));
    PCCMOC(f,19) = table(sum((Ch3o4_M(:,2).*Ch3o4_mask(:,1))/sum(Ch3o4_M(:,2))));
else end

%%%%Currently non-functioning math for % of overlapping pixels. It was
%%%%incorrectly being used to determine M1 and M2 Manders Coefficients.
% PCCMOC(f,8) = table(sum(((Ch1o2_M(:,1)>0).*(Ch1o2_M(:,2)>0)))/sum(Ch1o2_M(:,1)>0));
% PCCMOC(f,14) = table(sum(((Ch1o2_M(:,2)>0).*(Ch1o2_M(:,1)>0)))/sum(Ch1o2_M(:,2)>0));
% if Channels>2,
%     PCCMOC(f,9) = table(sum((Ch1o3_M(:,1).*(Ch1o3_M(:,2)>0)))/sum(Ch1o3_M(:,1)));
%     PCCMOC(f,15) = table(sum((Ch1o3_M(:,2).*(Ch1o3_M(:,1)>0)))/sum(Ch1o3_M(:,2)));
%     PCCMOC(f,11) = table(sum((Ch2o3_M(:,1).*(Ch2o3_M(:,2)>0)))/sum(Ch2o3_M(:,1)));
%     PCCMOC(f,17) = table(sum((Ch2o3_M(:,2).*(Ch2o3_M(:,1)>0)))/sum(Ch2o3_M(:,2)));
% else end
% if Channels>3,
%     PCCMOC(f,10) = table(sum((Ch1o4_M(:,1).*(Ch1o4_M(:,2)>0)))/sum(Ch1o4_M(:,1)));
%     PCCMOC(f,16) = table(sum((Ch1o4_M(:,2).*(Ch1o4_M(:,1)>0)))/sum(Ch1o4_M(:,2)));
%     PCCMOC(f,12) = table(sum((Ch2o4_M(:,1).*(Ch2o4_M(:,2)>0)))/sum(Ch2o4_M(:,1)));
%     PCCMOC(f,18) = table(sum((Ch2o4_M(:,2).*(Ch2o4_M(:,1)>0)))/sum(Ch2o4_M(:,2)));
%     PCCMOC(f,13) = table(sum((Ch3o4_M(:,1).*(Ch3o4_M(:,2)>0)))/sum(Ch3o4_M(:,1)));
%     PCCMOC(f,19) = table(sum((Ch3o4_M(:,2).*(Ch3o4_M(:,1)>0)))/sum(Ch3o4_M(:,2)));
% else end
%%%%%%%%%%%%%%%%%%%%%%%
   
%%%%This calculates and writes the PCC of each combination.
if SAVE_PCC > 0, cd(ImageDir);
    writematrix(Ch1o2_P,'Ch1o2_PCC.txt');
    if Channels>2, writematrix(Ch1o3_P,'Ch1o3_PCC.txt'); writematrix(Ch2o3_P,'Ch2o3_PCC.txt'); else end;
    if Channels>3, writematrix(Ch1o4_P,'Ch1o4_PCC.txt'); writematrix(Ch2o4_P,'Ch2o4_PCC.txt'); writematrix(Ch3o4_P,'Ch3o4_PCC.txt'); else end;
else end;

if isempty(Ch1o2_P) == 0, PCCMOC(f,2) = {corr(Ch1o2_P(:,1),Ch1o2_P(:,2))};
else end

if Channels>2,
    if isempty(Ch1o3_P) == 0, PCCMOC(f,3) = {corr(Ch1o3_P(:,1),Ch1o3_P(:,2))};
    else end
    if isempty(Ch2o3_P) == 0, PCCMOC(f,5) = {corr(Ch2o3_P(:,1),Ch2o3_P(:,2))};
    else end
else end
if Channels>3,
    if isempty(Ch1o4_P) == 0, PCCMOC(f,4) = {corr(Ch1o4_P(:,1),Ch1o4_P(:,2))};
    else end
    if isempty(Ch2o4_P) == 0, PCCMOC(f,6) = {corr(Ch2o4_P(:,1),Ch2o4_P(:,2))};
    else end
    if isempty(Ch3o4_P) == 0, PCCMOC(f,7) = {corr(Ch3o4_P(:,1),Ch3o4_P(:,2))};
    else end
else end

cd(userpath);
end

cd(Folder);
writetable(PCCMOC,'PCCMOC Analysis.xlsx','WriteVariableNames',true);

toc