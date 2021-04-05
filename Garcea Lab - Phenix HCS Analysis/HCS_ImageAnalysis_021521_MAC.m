clear; clc; tic; cd(userpath);
%% Editables %%

Folder = 'FILEPATH/'; 
%Where are the images located? IMPORTANT: Use format 'FILEPATH/'. The apostrophes and ending slash are important.

FigShow = 0; %Do you want to show the Figure with segmentation overlay? (1=yes, 0=no)
FigAnnotate = 0; %Do you want to show the number designation for nuclei and cell objects on the Figure? (1=yes, 0=no)

FigSave = 0; %Do you want to save the Figure that is generated during analysis? (1=yes, 0=no)
InFocus_ImageSave = 0; %Do you want to save the in-focus image hyperstack (XYC)? (1=yes, 0=no)

Channels = 3; %How many fluorescent channels are in the image?

CH_DAPI = 1; %Which channel corresponds to DAPI signal?

CH_CellMask = 3; %Which channel corresponds to CellMask (CM) signal?
CM_LowSizeFilt = 1000; %What is the minimum area of cell objects that you want to include in your analysis? (pixels^2)

TF_Analysis = 0; %Will do TF translocalization analysis (1=yes, 0=no)
    CH_TF = 3; %Which channel corresponds to transcription factor (TF) signal?

aSMA_Analysis = 1; %Will do gradient-based aSMA analysis (1=yes, 0=no)
    CH_aSMA = 2; %Which channel corresponds to aSMA signal?
    aSMA_ImageSave = 0; %Do you want to save an image of the aSMA channel and gradient image? (1=yes, 0=no)

SplitClusterResults = 0; %Will split results for each well into clusters (>MinClusterSize) and non-clusters (<MinClusterSize). (1=yes, 0=no)
    MinClusterSize = 8; %This dictates how "clusters" are identified for result sorting. 
%%The number refers to how many nuclei must be detected for a "group" to be called a "cluster.

%% Analysis Pre-Reqs and Metadata %%

cd(Folder);
srcFiles = dir('*.tiff');

for u = 1:length(srcFiles)
    srcFiles_HiddenFilt(u,1) = srcFiles(u).name(1) == '.'; %Identifies the "hidden" files that Mac creates when transferring data between Windows and Macs.
end
srcFiles = srcFiles(~srcFiles_HiddenFilt); %Filters out these hidden files.

START = 1;
FINISH = length(srcFiles);
for u = 1:length(srcFiles)
    srcFiles_FOV(u,1) = str2num(srcFiles(u).name(8:9));
    srcFiles_WELL(u,1) = string(srcFiles(u).name(1:6));
end
Fields = max(srcFiles_FOV);
Wells = numel(srcFiles)/Fields;

if FigShow == 1, figure,
else end

%% Analysis %%
for f = START:FINISH
    time(f,1).ElapsedSeconds = toc;
    
%     try
    
clc
filename = strcat(Folder,srcFiles(f).name);
progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
disp(progress2);

if f == START,
        cd(Folder); mkdir('Analysis'); cd(Folder);
else end

if progress < 10,
    disp('Estimated time remaining will display after 10% of images are analyzed...');
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
    
Results(f).FileName = srcFiles(f).name;
Results(f).Well = srcFiles_WELL{f};
Results(f).FieldofView = srcFiles_FOV(f);

I = bfopen(filename);

ResY = size(I{1,1}{1,1},1);
ResX = size(I{1,1}{1,1},2);
ZPlanes = (length(I{1,1})/Channels);
Blank3D = zeros(ResY,ResX,ZPlanes);
Slices = Channels*ZPlanes; 

%% Parsing Focal Plane and Creating In-Focus Image %%

disp('Generating in-focus images...');
[Ch1,Ch2,Ch3,Ch4,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4,ZFocus] = InFocusImage(I,ZPlanes,ResY,ResX,Channels,f);

Results(f).ZFocus = ZFocus(f,1);

if ZPlanes > 3
    if Results(f).ZFocus == 1 || Results(f).ZFocus == ZPlanes
    warning('No in-focus plane found. Continuing to next image.');
    pause(2);
    continue
    else
    end
else
end

if InFocus_ImageSave == 1,
    disp('Saving in-focus images...');
    if f == START,
        cd(Folder); cd Analysis; mkdir('InFocusImages'); cd InFocusImages;
    else cd(Folder); cd Analysis; cd InFocusImages; 
    end
        
    InFocus_IMAGE = zeros(ResY,ResX,Channels, 'uint16');  
    InFocus_IMAGE(:,:,1) = Focus_Ch1;
    InFocus_IMAGE(:,:,2) = Focus_Ch2;
    if Channels > 2, InFocus_IMAGE(:,:,3) = Focus_Ch3; else end
    if Channels > 3, InFocus_IMAGE(:,:,4) = Focus_Ch4; else end
    InFocus_IMAGENAME = strcat(srcFiles(f).name(1:9),' In-Focus Image.tiff');
    bfsave(InFocus_IMAGE(:,:,:),InFocus_IMAGENAME);
else
end

%% Cell Mask Segmentation %%

disp('Segmenting Cell Mask signal...');
if CH_CellMask == 1, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch1,CM_LowSizeFilt,ResY,ResX);
elseif CH_CellMask == 2, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch2,CM_LowSizeFilt,ResY,ResX);
elseif CH_CellMask == 3, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch3,CM_LowSizeFilt,ResY,ResX);
elseif CH_CellMask == 4, [CM_IndCells,CMseg_cc,CMseg_props,CMseg_CentroidXY,CM_Watershed_BW2,CM_Watershed_Perim] = CellMaskSegmentation(Focus_Ch4,CM_LowSizeFilt,ResY,ResX); 
else end

%% Nuclear Segmentation %%

disp('Segmenting DAPI signal...');
if CH_DAPI == 1, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch1,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 2, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch2,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 3, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch3,ResY,CM_Watershed_BW2);
elseif CH_DAPI == 4, [DAPI_75Percentile,DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Focus_Ch4,ResY,CM_Watershed_BW2); 
else end

%% Cellular Analysis %%

disp('Quantifying image properties for each cell/nucleus...');
if CH_DAPI == 1, DAPI = Focus_Ch1;
elseif CH_DAPI == 2, DAPI = Focus_Ch2;
elseif CH_DAPI == 3, DAPI = Focus_Ch3;
elseif CH_DAPI == 4, DAPI = Focus_Ch4;
else end

clearvars NearestNuc;
[Results_CellAnalysis,NearestNucDistanceFiltered] = CellularAnalysis(DAPI_75Percentile,DAPI_Watershed_BW2,Channels,CMseg_props,CM_IndCells,DAPI,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4);

clearvars NucleiNumber AdjustedNucleiNumber;
for q = 1:size(Results_CellAnalysis,2)
    NucleiNumber(q,1) = Results_CellAnalysis(q).NucleiNumber; 
    AdjustedNucleiNumber(q,1) = Results_CellAnalysis(q).AdjustedNucleiNumber;
end

%% TF Nuclear Translocalization Analysis %%

clearvars NucCytoRatio_sorted
if TF_Analysis == 1
    disp('Quantifying TF localization in each nucleus...');
    if CH_TF == 1, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch1,Results_CellAnalysis,Blank3D);
    elseif CH_TF == 2, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch2,Results_CellAnalysis,Blank3D);
    elseif CH_TF == 3, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch3,Results_CellAnalysis,Blank3D);
    elseif CH_TF == 4, [TFAnalysis,NucCytoRatio_sorted] = NucTranslocation(Focus_Ch4,Results_CellAnalysis,Blank3D);
    else
    end
else
end

%% aSMA Activation Analysis (work in progress) %%

clearvars Grad_Stats
if aSMA_Analysis == 1
    disp('Quantifying aSMA gradient intensities for each cell...');
    if CH_aSMA == 1, aSMA = Focus_Ch1;
    elseif CH_aSMA == 2, aSMA = Focus_Ch2;
    elseif CH_aSMA == 3, aSMA = Focus_Ch3;
    elseif CH_aSMA == 4, aSMA = Focus_Ch4;
    else
    end
    
    [Grad_Stats,aSMAGradImageTotal] = aSMAActivation(aSMA,Results_CellAnalysis);
else
end

if aSMA_Analysis == 1 && aSMA_ImageSave == 1,
    disp('Saving aSMA Figure...');
    if f == START,
        cd(Folder); cd Analysis; mkdir('aSMAImages'); cd aSMAImages;
    else cd(Folder); cd Analysis; cd aSMAImages;
    end
    aSMA_IMAGE = zeros(ResY,ResX,2, 'uint16');
    if CH_aSMA == 1, aSMA_IMAGE(:,:,1) = Focus_Ch1
    elseif CH_aSMA == 2, aSMA_IMAGE(:,:,1) = Focus_Ch2;
    elseif CH_aSMA == 3, aSMA_IMAGE(:,:,1) = Focus_Ch3;
    elseif CH_aSMA == 4, aSMA_IMAGE(:,:,1) = Focus_Ch4; 
    else end    
    aSMA_IMAGE(:,:,2) = aSMAGradImageTotal;    
    aSMA_IMAGENAME = strcat(srcFiles(f).name(1:9),' aSMA Gradient Image.tiff');
    bfsave(aSMA_IMAGE(:,:,:),aSMA_IMAGENAME);
else end
   
%% Figure %%

if FigShow > 0,
    disp('Generating Figure...');
totalfilteredcellperims = false(ResY,ResX);
for v = 1:size(Results_CellAnalysis,2)
    if v == 1, totalfilteredcellperims = Results_CellAnalysis(v).LogPerim;
    else totalfilteredcellperims = or(totalfilteredcellperims,Results_CellAnalysis(v).LogPerim);
    end
end

figimage1 = imoverlay(DAPI.*75,imdilate(DAPI_Watershed_Perim,strel('disk',1)),[0 0 1]);
figimage1 = imoverlay(figimage1,imdilate(totalfilteredcellperims,strel('disk',1)),[0.5 0 0]);
figimage2 = imoverlay(CellMask.*25,imdilate(totalfilteredcellperims,strel('disk',1)),[0.5 0 0]);
figimage2 = imoverlay(figimage2,imdilate(DAPI_Watershed_Perim,strel('disk',1)),[0 0 0.5]);
figimage3 = zeros(ResY,ResX,3);

C = [figimage1 figimage2 figimage3];
imshow(C); title('Segmentation and Analysis Summary');

dtclock = fix(clock);
scriptname = mfilename('fullpath');
scriptname_slash = strfind(scriptname,'/');
scriptname_str = scriptname(scriptname_slash(end)+1:end);
scriptname_underscore = strfind(scriptname_str,'_');
scriptname_str(scriptname_underscore) = ' ';
title('Segmentation and Analysis Summary');
figtext1 = ['Script: ' scriptname_str];
if dtclock(5)>10, figtext2 = ['Date and Time of Analysis: ' num2str(dtclock(2)) '/' num2str(dtclock(3)) '/' num2str(dtclock(1)) ' at ' num2str(dtclock(4)) ':' num2str(dtclock(5))];
else figtext2 = ['Date and Time of Analysis: ' num2str(dtclock(2)) '/' num2str(dtclock(3)) '/' num2str(dtclock(1)) ' at ' num2str(dtclock(4)) ':0' num2str(dtclock(5))];
end
figtext3 = ['Image Filename: ' srcFiles(f).name];
figtext4 = [' '];
figtext5 = ['In-Focus Z-Plane: ' num2str(Results(f).ZFocus)];
figtext6 = ['Nuclei Detected: ' num2str(sum(NucleiNumber,'all'))];
figtext7 = ['Nuclei Detected (Adjusted): ' num2str(sum(AdjustedNucleiNumber,'all'))];
figtext8 = ['Cell Objects Detected: ' num2str(size(Results_CellAnalysis,2))];
figtext_all = {figtext1,figtext2,figtext3,figtext4,figtext5,figtext6,figtext7,figtext8};
t = text(3000,500,figtext_all);
t.Color = [1 1 1];
t.FontSize = 12;

if FigAnnotate >0,
    for b = 1:size(Results_CellAnalysis,2)
        cmtext = text((ResX+Results_CellAnalysis(b).CMCentroid(1)),(Results_CellAnalysis(b).CMCentroid(2)),num2str(b));
        cmtext.Color = [1 1 0];
        cmtext.FontSize = 10;
        cmtext.FontWeight = 'bold';
        for n = 1:numel(Results_CellAnalysis(b).NuclearProps)
            nuctext = text(Results_CellAnalysis(b).NuclearProps(n).Centroid(1),Results_CellAnalysis(b).NuclearProps(n).Centroid(2),num2str(n));
            nuctext.Color = [1 1 0];
            nuctext.FontSize = 10;
            nuctext.FontWeight = 'bold';
        end
    end       
else
end

drawnow; hold off;
else
end

if FigSave == 1
    disp('Saving Figure...');
    if f == START,
        cd(Folder); cd Analysis; mkdir('FigureImages'); cd FigureImages;
    else cd(Folder); cd Analysis; cd FigureImages;
    end
    ax = gca;
    Fig_Name = append(filename(end-17:end-9),' Segmentation.tif');
    exportgraphics(ax,Fig_Name,'ContentType','image','Resolution','400');
else
% elseif f > START
%     cd(Folder); cd Analysis; cd FigureImages;
end

%% Results %%

clearvars CellularMeans CellularSums CMAreas aSMAGradientMean;
disp('Collating results...')

for r = 1:size(Results_CellAnalysis,2)
    CellularMeans(r,1) = Results_CellAnalysis(r).MeanCh1;
    CellularMeans(r,2) = Results_CellAnalysis(r).MeanCh2;
    if Channels>2, CellularMeans(r,3) = Results_CellAnalysis(r).MeanCh3; else end;
    if Channels>3, CellularMeans(r,4) = Results_CellAnalysis(r).MeanCh4; else end;
    CellularSums(r,1) = Results_CellAnalysis(r).SumCh1;
    CellularSums(r,2) = Results_CellAnalysis(r).SumCh2;
    if Channels >2, CellularSums(r,3) = Results_CellAnalysis(r).SumCh3; else end;
    if Channels >3, CellularSums(r,4) = Results_CellAnalysis(r).SumCh4; else end;
    CMAreas(r,1) = Results_CellAnalysis(r).CellArea;    
       
    if aSMA_Analysis == 1, aSMAGradientMean(r,1) = Grad_Stats(r).MeanGradientValue;
    else
    end
end

if ZPlanes > 3 && Results(f).ZFocus == 1 || Results(f).ZFocus == ZPlanes
    Results(f).NucNearestNeighborDistance = "";
    Results(f).TotalNuclei = 0;
    Results(f).NucleiPerGroup = 0;
    Results(f).AdjustedTotalNuclei = 0;
    Results(f).AdjustedNucGroupSizes = "";
    Results(f).MeanIntensities = ""; %Each column is a channel.
    Results(f).SumIntensities = ""; %Each column is a channel.
    Results(f).CMNumber = 0;
    Results(f).CMAreas = "";
    Results(f).CMTotalArea = 0;
    if TF_Analysis==1, Results(f).TFNucCytoRatios = 0;
    else end
    if aSMA_Analysis>0, Results(f).aSMAGradientMeans = 0;
    else end
else
    Results(f).NucNearestNeighborDistance = NearestNucDistanceFiltered;
    Results(f).TotalNuclei = sum(NucleiNumber,'all');
    Results(f).NucleiPerGroup = NucleiNumber;
    Results(f).AdjustedTotalNuclei = sum(AdjustedNucleiNumber,'all');
    Results(f).AdjustedNucGroupSizes = AdjustedNucleiNumber;
    Results(f).MeanIntensities = CellularMeans; %Each column is a channel.
    Results(f).SumIntensities = CellularSums; %Each column is a channel.
    Results(f).CMNumber = size(Results_CellAnalysis,2);
    Results(f).CMAreas = CMAreas;
    Results(f).CMTotalArea = sum(CMAreas);
    if TF_Analysis>0, Results(f).TFNucCytoRatios = NucCytoRatio_sorted;
    else end
    if aSMA_Analysis>0, Results(f).aSMAGradientMeans = aSMAGradientMean;
    else end
end

%     catch
%          warning('An error occurred during analysis. Saving Results and skipping to next image.'); pause(2);
%          Results = Results(START:f);
%          cd(Folder); cd Analysis;
%          save('AnalysisResults.mat','Results', '-v7.3');
%     end
end %End of Analysis Loop

cd(Folder); cd Analysis;
save('AnalysisResults.mat','Results', '-v7.3');

    %% Result Sorting By Well%%

disp('Finished analyzing images, now sorting results by well...');
[Results_Sorted,WellNumber] = WellSort(Channels,Results,ZPlanes,TF_Analysis,aSMA_Analysis);
cd(Folder); cd Analysis;
save('AnalysisResultsSortedByWell.mat','Results_Sorted','-v7.3');
    
    %% Sorting Results by Cluster Status %%

if SplitClusterResults == 1
    disp('Sorting results by cluster status...');
    [Results_Sorted_InClusters,Results_Sorted_OutClusters] = ClusterSort(Channels,Results_Sorted,WellNumber,MinClusterSize,TF_Analysis,aSMA_Analysis);

cd(Folder); cd Analysis;
    save('AnalysisResultsSortedByWellandWithinClusters.mat','Results_Sorted_InClusters','-v7.3');
    save('AnalysisResultsSortedByWellandOutsideClusters.mat','Results_Sorted_OutClusters','-v7.3');
else end

toc