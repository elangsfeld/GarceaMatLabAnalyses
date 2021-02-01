clear all; clc; tic; cd(userpath);

%% Editables %%

Folder = 'F:\DATA\08.02.19 EdU P-C with RA, NG18, and 808A\';

Channels = 3; %How many fluorescent channels are present in the images? (Values must = 2-4).
ShowFig = 1; %Option to display the Figure of the image and VRC IDs (1=yes, 0=no).
SaveFig = 1; %Option to save the Figure image as a PDF (1=yes,0=no).
ImageSave_InFocus = 0; %Option to save the "in-focus" plane of each channel in a hyperstack (XYC) (1=yes,0=no).
ImageSave_PostThreshold = 0; %Option to save the post-threshold images of each channel in a hyperstack (XYZC) (1=yes,0=no).

Percentile_Thresholds.Lower(1:4) = [99.75 99.75 99.75 99.75]; %What is the lower bound percentile-based threshold you want to use for analysis?
Percentile_Thresholds.Upper(1:4) = [100 100 100 100]; %What is the upper bound percentile-based threshold you want to use for analysis?

DBSCAN_Rad = [15 15 15 15];
DBSCAN_MinPts = [50 50 50 50];

VRC_ID = [1 0 0 0]; %Chooses which channel(s) to use for VRC identification (can be up to 2).
RS_Rad = [40 40 40 40]; 
VRC_RSseed = 1; %Chooses which channel to use for VRC point seeding (for filtering rangesearch).

Analysis_Volume = 0;
Analysis_PCC = 0;

%% Pre-Analysis Things %%

cd(Folder);
srcFiles = dir('*.nd2');

START = 1;
FINISH = numel(srcFiles);

if sum(VRC_ID(:)) < 2,
    VRC_Idx2(1:4,2) = [NaN;NaN;NaN;NaN];
else end

%% Pre-Analysis Things for Each Image %%

for f = START:FINISH
    time(f,1).ElapsedSeconds = toc;
    
    clc
    filename = strcat(Folder,srcFiles(f).name);
    progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2)
       
    if f == START,
        cd(Folder); mkdir('AnalysisResults'); cd(Folder);
    else end
    
    if progress < 10,
        disp('Estimated time remaining will display after 10% of images have been analyzed...');
    else
        time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/((FINISH-START+1)-(FINISH-f));
        time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH-START+1);
        time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
        time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
        time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
        estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
        estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
        disp(estimate);
        disp(estimate2);
    end
    
    Results(f).FileName = srcFiles(f).name(1:end-4);
    
    RawImage = bfopen(filename);
    
    Res = length(RawImage{1,1}{1,1});
    Slices = (length(RawImage{1,1})/Channels);
    Blank3D = zeros(Res,Res,Slices);

%% Separate Channels and Parse Focal Plan %%
disp('Splitting Channels and Creating In-Focus Images...');
if f==START, FocusInfo = struct(); else end
[I,FocusInfo] = SplitND2Stack(Channels,Res,RawImage,f,FocusInfo);

Results(f).ZFocus = mode(FocusInfo(f).MaxPlane);

if Results(f).ZFocus == 1 || Results(f).ZFocus == Slices
    warning('No in-focus plane found. Continuing to next image.');
    pause(2);
continue
else
end

if ImageSave_InFocus == 1,
    disp('Saving in-focus images...');
    if f == START
        cd(Folder); cd AnalysisResults; mkdir('InFocusImages'); cd InFocusImages;
    else cd(Folder); cd AnalysisResults; cd InFocusImages;
    end
    IMAGE_InFocus = zeros(Res,Res,Channels,'uint16');
    IMAGE_InFocus(:,:,1) = I.Ch1(:,:,Results(f).ZFocus);
    IMAGE_InFocus(:,:,2) = I.Ch2(:,:,Results(f).ZFocus);
    if Channels > 2, IMAGE_InFocus(:,:,3) = I.Ch3(:,:,Results(f).ZFocus); else end;
    if Channels > 3, IMAGE_InFocus(:,:,4) = I.Ch4(:,:,Results(f).ZFocus); else end;
    ImageName_InFocus = strcat(Results(f).FileName,' In-Focus Image.tiff');
    bfsave(IMAGE_InFocus(:,:,:),ImageName_InFocus);
else
end

%% Percentile-Based Thresholding %%
disp('Calculating percentile-based thresholds and applying to images...');
[Percentile_Values,BWMask,Seg] = PercentileThresholding(f,I,Channels,Percentile_Thresholds);

if ImageSave_PostThreshold == 1
    disp('Saving post-threshold images...');
    if f == START
        cd(Folder); cd AnalysisResults; mkdir('PostThresholdImages'); cd PostThresholdImages;
    else cd(Folder); cd AnalysisResults; cd PostThresholdImages;
    end

    ImageName_PostThreshCh1 = strcat(Results(f).FileName,' Post-Threshold Image Ch1.tiff');
    ImageName_PostThreshCh2 = strcat(Results(f).FileName,' Post-Threshold Image Ch2.tiff');
    bfsave(Seg.Ch1,ImageName_PostThreshCh1);
    bfsave(Seg.Ch2,ImageName_PostThreshCh2);
    if Channels > 2
        ImageName_PostThreshCh3 = strcat(Results(f).FileName,' Post-Threshold Image Ch3.tiff');
        bfsave(Seg.Ch3,ImageName_PostThreshCh3);
    else end
    if Channels > 3
        ImageName_PostThreshCh4 = strcat(Results(f).FileName,' Post-Threshold Image Ch4.tiff');
        bfsave(Seg.Ch4,ImageName_PostThreshCh4);
    else end
else
end

%% DBSCAN analysis from binary masks %%
disp('Performing DBSCAN analysis to identify signal clusters...');
if f==START, ClusterInfo = struct(); DBSCAN = struct(); else end
[ClusterInfo,DBSCAN] = VRC_DBSCAN(f,DBSCAN_Rad,DBSCAN_MinPts,Channels,VRC_ID,FocusInfo,BWMask,ClusterInfo,DBSCAN);

%% RangeSearch Nearest Neighbor to combine VRCs IDed by different channels, and create expanded Bounding Boxes %%
disp('Performing RangeSearch analysis to identify nearest-neighbor VRC(s) for subsequent analyses...');
if f==START, RSInfo = []; else end
[RSInfo] = VRC_RS(f,RS_Rad,VRC_ID,VRC_RSseed,ClusterInfo,RSInfo);
close all;

%% Draw VRCs on image and save nuclear image as PDF %%

for slash = 1:length(strfind(filename,'\'))
    if slash ==1, filename2 = extractAfter(filename,'\');
    else filename2 = extractAfter(filename2,'\');
    end
end
filename2 = filename2(1:end-4); 
filename3 = append(filename2,' Figure Image.pdf');
Folder_Outputs = append(Folder,'Outputs\');
    
if f==START
    cd(Folder); cd AnalysisResults; mkdir('FigureImages'); cd FigureImages;
else cd(Folder); cd AnalysisResults; cd FigureImages;
end
    
if ShowFig == 1,
    figure,
    if VRC_ID == [1,0,0,0], imshow(Seg.Ch1(:,:,FocusInfo(f).MaxPlane(1,1)).*3);
    elseif VRC_ID == [0,1,0,0], imshow(Seg.Ch2(:,:,FocusInfo(f).MaxPlane(1,2)).*3);
    elseif VRC_ID == [0,0,1,0], imshow(Seg.Ch3(:,:,FocusInfo(f).MaxPlane(1,3)).*3);
    elseif VRC_ID == [0,0,0,1], imshow(Seg.Ch4(:,:,FocusInfo(f).MaxPlane(1,4)).*3);
    elseif VRC_ID == [1,1,0,0], imshow(imfuse(Seg.Ch1(:,:,FocusInfo(f).MaxPlane(1,1)).*3,Seg.Ch2(:,:,FocusInfo(f).MaxPlane(1,2)).*3,'ColorChannels',[2 1 1]));
    elseif VRC_ID == [1,0,1,0], imshow(imfuse(Seg.Ch1(:,:,FocusInfo(f).MaxPlane(1,1)).*3,Seg.Ch3(:,:,FocusInfo(f).MaxPlane(1,3)).*3,'ColorChannels',[1 2 0]));
    elseif VRC_ID == [1,0,0,1], imshow(imfuse(Seg.Ch1(:,:,FocusInfo(f).MaxPlane(1,1)).*3,Seg.Ch4(:,:,FocusInfo(f).MaxPlane(1,4)).*3,'ColorChannels',[2 1 1]));
    elseif VRC_ID == [0,1,1,0], imshow(imfuse(Seg.Ch2(:,:,FocusInfo(f).MaxPlane(1,2)).*3,Seg.Ch3(:,:,FocusInfo(f).MaxPlane(1,3)).*3,'ColorChannels',[2 1 0]));
    elseif VRC_ID == [0,1,0,1], imshow(imfuse(Seg.Ch2(:,:,FocusInfo(f).MaxPlane(1,2)).*3,Seg.Ch4(:,:,FocusInfo(f).MaxPlane(1,4)).*3,'ColorChannels',[1 2 1]));
    elseif VRC_ID == [0,0,1,1], imshow(imfuse(Seg.Ch3(:,:,FocusInfo(f).MaxPlane(1,3)).*3,Seg.Ch4(:,:,FocusInfo(f).MaxPlane(1,4)).*3,'ColorChannels',[1 2 1]));
    else end
    
    hold on
    
    for a = 1:size(RSInfo(f).ExpandedBB,1)
        rectangle('Position',[RSInfo(f).ExpandedBB(a).MinX RSInfo(f).ExpandedBB(a).MinY (RSInfo(f).ExpandedBB(a).MaxX-RSInfo(f).ExpandedBB(a).MinX) (RSInfo(f).ExpandedBB(a).MaxY-RSInfo(f).ExpandedBB(a).MinY)],'EdgeColor','y'); %[x y w h]
    end
    
    hold off
   
    if SaveFig == 1,
    ax = gca;
    exportgraphics(ax,filename3,'ContentType','vector','Resolution','500');
    else
    end

else
end

% %% Dissect individual VRCs and save images (if drawVRCimages == 1).
% 
%     cd(Folder_Outputs);
%     Folder_VRCs=append(filename2,' VRCs\');
% 
% if ShowFig == 1,
%     for t = 1:size(VRC_crops,1),
%         if t == 1, mkdir(Folder_VRCs); cd(Folder_VRCs);
%             delete *.tif;
%         else end
%         
%         VRC_subname = num2str(t);
%         if t <10, VRC_filename = append(filename2,' VRC 0',VRC_subname,'.tiff');
%         else VRC_filename = append(filename2,' VRC ',VRC_subname,'.tiff');
%         end
%         
%         VRCs{t,1}(:,:,:,1) = Ch1_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
%         VRCs{t,1}(:,:,:,2) = Ch2_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
%         if Channels >2, VRCs{t,1}(:,:,:,3) = Ch3_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
%         if Channels >3, VRCs{t,1}(:,:,:,4) = Ch4_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
%                 
%         bfsave(VRCs{t,1}(:,:,:,:),VRC_filename,'dimensionOrder','XYZTC');
%     end
%     
% else 
%     for t = 1:size(VRC_crops,1),       
%         VRCs{t,1}(:,:,:,1) = Ch1_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
%         VRCs{t,1}(:,:,:,2) = Ch2_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
%         if Channels >2, VRCs{t,1}(:,:,:,3) = Ch3_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
%         if Channels >3, VRCs{t,1}(:,:,:,4) = Ch4_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
%     end
%  
% end
% 
% %% Volume Analysis for each VRC %%
% % VRCs{#,1}(X,Y,Z,C).
% 
% if Analysis_Volume == 1
%     
% clearvars VRC_Volumes;
% VRC_Volumes = table('Size',[size(VRCs,1) length(varNames_Volumes)],'VariableTypes',varTypes_Volumes,'VariableNames',varNames_Volumes);
% 
% cd(Folder_Outputs);
% if f == START, mkdir('Volume Analysis'); cd('Volume Analysis');
%     else cd('Volume Analysis'); 
% end
% 
% for h = 1:size(VRCs,1),
%     
%     Volumes_filename = append(filename2,' Volumes.xlsx');
% 
%     for u = 1:size(VRCs,1),
%         VRC_Volumes.Ch1Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,1)>0)));
%         VRC_Volumes.Ch2Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,2)>0)));
%         if Channels>2, VRC_Volumes.Ch3Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,3)>0))); else, VRC_Volumes.Ch3Pix(:) = 0; end
%         if Channels>3, VRC_Volumes.Ch4Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,4)>0))); else, VRC_Volumes.Ch4Pix(:) = 0; end
%     end
% end
% 
% writetable(VRC_Volumes,Volumes_filename,'WriteVariableNames',true); cd(Folder);
% 
% else end
% %% PCC Analysis for each VRC %%
% 
% if Analysis_PCC ==1,
% 
% clearvars VRC_PCCs;
% VRC_PCCs = table('Size',[size(VRCs,1) length(varNames_PCCs)],'VariableTypes',varTypes_PCCs,'VariableNames',varNames_PCCs);
% 
% cd(Folder_Outputs);
% if f == START, mkdir('PCC Analysis'); cd('PCC Analysis');
%     else cd('PCC Analysis');
% end
% 
% for hh = 1:size(VRCs,1),
%     
%     PCCs_filename = append(filename2,' PCCs.xlsx');
%     VRC_area = (VRC_crops.VRC_MaxY(1) - VRC_crops.VRC_MinY(1))*((VRC_crops.MaxX(1) - VRC_crops.VRC_minX(1)));
%     
%     for uu = 1:size(VRCs,1),
%         Ch1_re = reshape(Ch1_seg,[VRC_area 1]);
%         
%     end
% end
% else end
end
close all %Closes all open figures that are made from "gscatter". Not sure if necessarily the best way, but it works.
toc




