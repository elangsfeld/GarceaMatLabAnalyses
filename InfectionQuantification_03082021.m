
clear all; clc; tic; cd(userpath);

%% Editables %%

Folder = 'D:\6515\';
%Format: 'FILEPATH\'

FigShow = 1; %Do you want to display an image of nuclear segmentation and TAg signal? (1=yes,0=no).
    FigSave = 0; %Do you want to save a .jpeg image of the figure? (1=yes,0=no).

Channels = 3; %How many channels/wavelengths are in the data set?

Ch_DAPI = 1; %Which channel number corresponds to DAPI/Hoechst signal?
    LowAreaBound = 10; %What is the smallest DAPI object you want to detect? This value could change if camera binning values change.
    HighAreaBound = 500; %What is the largest DAPI object you want to detect? This value could change if camera binning values change.

Analysis_TAg = 1; %Do you want to analyze TAg intensity? (1=yes,0=no).
    Ch_TAg = 2; %Which channel number corresponds to TAg signal?
    Threshold_TAg = 2000; %Where do you want the threshold for infected and uninfected cells to be called, based on TAg signal?
    FigShow_TAg = 1; %Do you want to display the TAg image in the figure? (1=yes,0=no).

Analysis_VP1 = 1; %Do you want to analyze VP1 intensity? (1=yes,0=no).
    Ch_VP1 = 3; %Which channel number corresponds to VP1 signal?
    Threshold_VP1 = 2000; %Where do you want the threshold for infected and uninfected cells to be called, based on VP1 signal?
    FigShow_VP1 = 1; %Do you want to display the VP1 image in the figure? (1=yes,0=no).
  

%% File Sorting and Property Collection %%
disp('Collecting file information and sorting image channels...');

cd(Folder);
srcFiles = dir('*.tif');

START = 1;
FINISH = length(srcFiles)/Channels;

for i = START:length(srcFiles)
    if numel(strfind(srcFiles(i).name,'w1')) > 0
        SortIdx(i,1) = 1;
    elseif numel(strfind(srcFiles(i).name,'w2')) > 0
        SortIdx(i,2) = 1;
    elseif Channels>2 && numel(strfind(srcFiles(i).name,'w3'))>0
        SortIdx(i,3) = 1;
    elseif Channels>3 && numel(strfind(srcFiles(i).name,'w4'))>0
    else
    end
end

srcFilesSorted.Ch1 = srcFiles(logical(SortIdx(:,1)));
srcFilesSorted.Ch2 = srcFiles(logical(SortIdx(:,2)));
if Channels>2, srcFilesSorted.Ch3 = srcFiles(logical(SortIdx(:,3))); else end
if Channels>3, srcFilesSorted.Ch4 = srcFiles(logical(SortIdx(:,4))); else end

for f = 1:size(srcFilesSorted.Ch1,1) 
    FileInfo.Ch1(f,1).Filename = srcFilesSorted.Ch1(f,1).name; 
    namesort(f,1).underscoreIdx = strfind(FileInfo.Ch1(f,1).Filename,'_');
    namesort(f,1).well = namesort(f,1).underscoreIdx(1,1)+1;
    namesort(f,1).site = namesort(f,1).underscoreIdx(1,2)+2;
    namesort(f,1).channel = namesort(f,1).underscoreIdx(1,3)+2;
    FileInfo.Ch1(f,1).Well = FileInfo.Ch1(f,1).Filename(namesort(f,1).well:(namesort(f,1).well+2));
    FileInfo.Ch1(f,1).Site = FileInfo.Ch1(f,1).Filename(namesort(f,1).site:namesort(f,1).site);
    FileInfo.Ch1(f,1).Channel = FileInfo.Ch1(f,1).Filename(namesort(f,1).channel:namesort(f,1).channel); 
    
    FileInfo.Ch2(f,1).Filename = srcFilesSorted.Ch2(f,1).name;
    FileInfo.Ch2(f,1).Well = FileInfo.Ch2(f,1).Filename(namesort(f,1).well:(namesort(f,1).well+2));
    FileInfo.Ch2(f,1).Site = FileInfo.Ch2(f,1).Filename(namesort(f,1).site:namesort(f,1).site);
    FileInfo.Ch2(f,1).Channel = FileInfo.Ch2(f,1).Filename(namesort(f,1).channel:namesort(f,1).channel);
    
    if Channels>2
        FileInfo.Ch3(f,1).Filename = srcFilesSorted.Ch3(f,1).name;
        FileInfo.Ch3(f,1).Well = FileInfo.Ch3(f,1).Filename(namesort(f,1).well:(namesort(f,1).well+2));
        FileInfo.Ch3(f,1).Site = FileInfo.Ch3(f,1).Filename(namesort(f,1).site:namesort(f,1).site);
        FileInfo.Ch3(f,1).Channel = FileInfo.Ch3(f,1).Filename(namesort(f,1).channel:namesort(f,1).channel);
    else end
    
    if Channels>3
        FileInfo.Ch4(f,1).Filename = srcFilesSorted.Ch4(f,1).name;
        FileInfo.Ch4(f,1).Well = FileInfo.Ch4(f,1).Filename(namesort(f,1).well:(namesort(f,1).well+2));
        FileInfo.Ch4(f,1).Site = FileInfo.Ch4(f,1).Filename(namesort(f,1).site:namesort(f,1).site);
        FileInfo.Ch4(f,1).Channel = FileInfo.Ch4(f,1).Filename(namesort(f,1).channel:namesort(f,1).channel);
    else end
    
    SiteNum(f,1) = str2num(FileInfo.Ch1(f,1).Site); 
end
    SiteMax = max(SiteNum);
  
if FigShow == 1, figure; else end
for f = START:FINISH
    time(f,1).ElapsedSeconds = toc;
  

    %    try
    clc
    progress = (((FINISH-START+1)-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2);
    
    if f == START
        cd(Folder); mkdir('Analysis'); cd(Folder);
    else
        cd(Folder);
    end

    if progress < 10
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
       
    I.Ch1 = imread(FileInfo.Ch1(f,1).Filename);
    I.Ch2 = imread(FileInfo.Ch2(f,1).Filename);
    if Channels > 2, I.Ch3 = imread(FileInfo.Ch3(f,1).Filename); else end
    if Channels > 3, I.Ch4 = imread(FileInfo.Ch4(f,1).Filename); else end
    
    Res = size(I.Ch1,1);
        
    Results(f,1).Filename = FileInfo.Ch1(f).Filename;
    Results(f,1).Well = FileInfo.Ch1(f).Well;
    Results(f,1).Site = FileInfo.Ch1(f).Site;
    
    %% Nuclei Segmentation %%
    disp('Segmenting nuclei based on DAPI signal...');    
    if Ch_DAPI == 1, DAPI.Image = I.Ch1;
    elseif Ch_DAPI == 2, DAPI.Image = I.Ch2;
    elseif Ch_DAPI == 3, DAPI.Image = I.Ch3;
    elseif Ch_DAPI == 4, DAPI.Image = I.Ch4;
    else end
    
    DAPI.bg = imopen(DAPI.Image,strel('disk',round(Res/25)));
    DAPI.Image2 = DAPI.Image - DAPI.bg;
    DAPI.BW = imbinarize(DAPI.Image2);
    DAPI.BW2 = imfill(DAPI.BW,'holes');
    DAPI.BW2 = xor(bwareaopen(DAPI.BW2,LowAreaBound),bwareaopen(DAPI.BW2,HighAreaBound));
    
    DAPI.cc = bwconncomp(DAPI.BW,4);
    
    figimage1 = imoverlay(imadjust(DAPI.Image2),imdilate(bwperim(DAPI.BW2),strel('disk',1)),[1 0 0]);
    
    RP(f,1).Ch1 = regionprops(DAPI.cc,I.Ch1,'Area','MeanIntensity');
    RP(f,1).Ch2 = regionprops(DAPI.cc,I.Ch2,'Area','MeanIntensity');
    if Channels>2, RP(f,1).Ch3 = regionprops(DAPI.cc,I.Ch3,'Area','MeanIntensity'); else end
    if Channels>3, RP(f,1).Ch4 = regionprops(DAPI.cc,I.Ch4,'Area','MeanIntensity'); else end
    
    %% TAg Signal Analysis %%
    if Analysis_TAg == 1
        disp('Analyzing TAg signal intensity within each nucleus...');

        if f>START, clearvars TAg; else end

        if Ch_TAg == 1
            TAg.Image = I.Ch1;
            TAg.RP = RP(f,1).Ch1;
        elseif Ch_TAg == 2
            TAg.Image = I.Ch2;
            TAg.RP = RP(f,1).Ch2;
        elseif Ch_TAg == 3 
            TAg.Image = I.Ch3;
            TAg.RP = RP(f,1).Ch3;
        elseif Ch_TAg == 4
            TAg.Image = I.Ch4;
            TAg.RP = RP(f,1).Ch4;
        else
        end

        for t = 1:length(TAg.RP)
            TAg.MeanIntensity(t,1) = TAg.RP(t,1).MeanIntensity;
            TAg.AboveThresh(t,1) = TAg.MeanIntensity(t,1) > Threshold_TAg;
        end
        TAg.NumberAboveThresh = sum(TAg.AboveThresh,'all');
        TAg.NumberTotal = length(TAg.RP);
        TAg.PercentInfected = TAg.NumberAboveThresh/TAg.NumberTotal*100;

        if FigShow == 1 && FigShow_TAg == 1
            figimage2 = imoverlay(imadjust(TAg.Image),bwperim(DAPI.BW2),[1 0 0]);
        else
        end
    else
    end
    
    %% VP1 Signal Analysis %%
    if Analysis_VP1 == 1
        disp('Analyzing VP1 signal intensity within each nucleus...');

        if f>START, clearvars VP1; else end

        if Ch_VP1 == 1
            VP1.Image = I.Ch1;
            VP1.RP = RP(f,1).Ch1;
        elseif Ch_VP1 == 2
            VP1.Image = I.Ch2;
            VP1.RP = RP(f,1).Ch2;
        elseif Ch_VP1 == 3
            VP1.Image = I.Ch3;
            VP1.RP = RP(f,1).Ch3;
        elseif Ch_VP1 == 4
            VP1.Image = I.Ch4;
            VP1.RP = RP(f,1).Ch4;
        else
        end

        for v = 1:length(VP1.RP)
            VP1.MeanIntensity(v,1) = VP1.RP(v,1).MeanIntensity;
            VP1.AboveThresh(v,1) = VP1.MeanIntensity(v,1) > Threshold_VP1;
        end
        VP1.NumberAboveThresh = sum(VP1.AboveThresh,'all');
        VP1.NumberTotal = length(VP1.RP);
        VP1.PercentInfected = VP1.NumberAboveThresh/VP1.NumberTotal*100;

        if FigShow ==1 && FigShow_VP1 == 1
            figimage3 = imoverlay(imadjust(VP1.Image),bwperim(DAPI.BW2),[1 0 0]);
        else
        end
    else
    end
    
    %% Figure Generation and Saving %%
    if FigShow == 1
        disp('Generating Figure...');
        if FigShow_TAg == 1 && FigShow_VP1 == 1
            C = [figimage1 figimage2 figimage3];
        elseif FigShow_TAg == 1 && FigShow_VP1 == 0
            C = [figimage1 figimage2];
        elseif FigShow_TAg == 0 && FigShow_VP1 == 1
            C = [figimage1 figimage3];
        elseif FigShow_TAg == 0 && FigShow_VP1 == 0
            C = figimage1;
        else
        end
    else
    end
    
    if FigShow == 1, imshow(C); else end
    
    if FigSave == 1
        disp('Saving Figure Image...');
        if f == START
            cd(Folder); cd Analysis; mkdir('AnalysisImages'); cd AnalysisImages;
        else cd(Folder); cd Analysis; cd AnalysisImages;
        end
        ax = gca;
        FigName = FileInfo.Ch1(f,1).Filename;
        FigName = append(FileInfo.Ch1(f,1).Filename(1:strfind(FigName,'_w')),'Segmentation Image.jpg');
        exportgraphics(ax,FigName,'ContentType','image','Resolution','400');
    else
    end
    
    %% Results Collection %%
    disp('Updating Results...');
    
    
    
    Results_All(f,1).Filename = FileInfo.Ch1(f,1).Filename(1:(strfind(FileInfo.Ch1(f,1).Filename,'_w')-1));
    Results_All(f,1).Well = FileInfo.Ch1(f,1).Well;
    Results_All(f,1).Site = str2num(FileInfo.Ch1(f,1).Site);
    Results_All(f,1).NumberNucleiDetected = DAPI.cc.NumObjects;
    if Analysis_TAg == 1
        Results_All(f,1).NumberTAgPositive = TAg.NumberAboveThresh;
        Results_All(f,1).PercentTAgPositive = TAg.PercentInfected;
    else end
    if Analysis_VP1 == 1
        Results_All(f,1).NumberVP1Positive = VP1.NumberAboveThresh;
        Results_All(f,1).PercentVP1Positive = VP1.PercentInfected;
    else end
    
end

%% Well/Plate Sorting %%
disp('Sorting results by well and saving .mat files...');
Idx_WellStart = [1:SiteMax:FINISH];
Idx_WellEnd = [SiteMax:SiteMax:FINISH];
WellNumber = length(Idx_WellStart);

for w = 1:WellNumber
    for f = Idx_WellStart(w):Idx_WellEnd(w)
        if f == Idx_WellStart(w)
            Results_WellSorted(w,1).Well = Results_All(f,1).Well;
            Results_WellSorted(w,1).NumberNucleiDetected = Results_All(f,1).NumberNucleiDetected;
            if Analysis_TAg == 1
                Results_WellSorted(w,1).NumberTAgPositive = Results_All(f,1).NumberTAgPositive;
            else end
            if Analysis_VP1 == 1
                Results_WellSorted(w,1).NumberVP1Positive = Results_All(f,1).NumberVP1Positive;
            else end
        else
            Results_WellSorted(w,1).NumberNucleiDetected = Results_WellSorted(w,1).NumberNucleiDetected + Results_All(f,1).NumberNucleiDetected;
            if Analysis_TAg ==1
                Results_WellSorted(w,1).NumberTAgPositive = Results_WellSorted(w,1).NumberTAgPositive + Results_All(f,1).NumberTAgPositive;
            else end
            if Analysis_VP1 == 1
                Results_WellSorted(w,1).NumberVP1Positive = Results_WellSorted(w,1).NumberVP1Positive + Results_All(f,1).NumberVP1Positive;
            else end
        end
    end
    
    if Analysis_TAg == 1, Results_WellSorted(w,1).PercentTAgPositive = Results_WellSorted(w,1).NumberTAgPositive/Results_WellSorted(w,1).NumberNucleiDetected*100;
    else end
    if Analysis_VP1 == 1, Results_WellSorted(w,1).PercentVP1Positive = Results_WellSorted(w,1).NumberVP1Positive/Results_WellSorted(w,1).NumberNucleiDetected*100;
    else end   
end

cd(Folder); cd Analysis; mkdir('InfectivityQuantificationResults'); cd InfectivityQuantificationResults;
save('Infectivity Analysis Results.mat','Results_All','-v7.3');
save('Infectivity Analysis Results, Sorted By Well.mat','Results_WellSorted','-v7.3');