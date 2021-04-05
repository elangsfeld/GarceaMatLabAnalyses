clear; clc; tic; cd(userpath);
%% Editables %%

Folder = 'F:\DATA\02.26.21 Phenix Imaging of iPOND Abs and others\Images\ImageStacks\'; 
%Where are the images located? IMPORTANT: Use format 'FILEPATH\'. The apostrophes and ending slash are important.

FigShow = 0; %Do you want to show the Figure with segmentation overlay? (1=yes, 0=no)
    FigAnnotate = 1; %Do you want to show the number designation for nuclei and cell objects on the Figure? (1=yes, 0=no)

FigSave = 0; %Do you want to save the Figure that is generated during analysis? (1=yes, 0=no)
InFocus_ImageSave = 0; %Do you want to save the in-focus image hyperstack (XYC)? (1=yes, 0=no)

Channels = 4; %How many fluorescent channels are in the image?

CH_DAPI = 1; %Which channel corresponds to DAPI signal?

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

if CH_DAPI == 1, DAPI = Focus_Ch1;
elseif CH_DAPI == 2, DAPI = Focus_Ch2;
elseif CH_DAPI == 3, DAPI = Focus_Ch3;
elseif CH_DAPI == 4, DAPI = Focus_Ch4;
else warning('Check value of CH_DAPI. Must equal 1-4...');
end

%% Nuclear Segmentation %%

disp('Segmenting DAPI signal...');
[DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(DAPI,ResY);

%% Nuclear Analysis %%

disp('Quantifying image properties for each nucleus...');

clearvars NearestNuc;
[Results_NuclearAnalysis] = NuclearAnalysis(DAPI_Watershed_BW2,Channels,DAPI,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4);
   
%% Figure %%

if FigShow > 0
    disp('Generating Figure...');

figimage1 = imoverlay(DAPI.*5,imdilate(DAPI_Watershed_Perim,strel('disk',1)),[0 0 1]);
figimage2 = zeros(ResY,ResX,3);

C = [figimage1 figimage2];
imshow(C); title('Segmentation and Analysis Summary');

dtclock = fix(clock);
scriptname = mfilename('fullpath');
scriptname_slash = strfind(scriptname,'\');
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
figtext5 = ['Nuclei Detected: ' num2str(sum(Results_NuclearAnalysis.NucleiNumber,'all'))];
figtext_all = {figtext1,figtext2,figtext3,figtext4,figtext5};
t = text(2500,300,figtext_all);
t.Color = [1 1 1];
t.FontSize = 16;

if FigAnnotate >0,
    for n = 1:numel(Results_NuclearAnalysis.NuclearProps)
        nuctext = text(Results_NuclearAnalysis.NuclearProps(n).Centroid(1),Results_NuclearAnalysis.NuclearProps(n).Centroid(2),num2str(n));
        nuctext.Color = [0.75 0 0];
        nuctext.FontSize = 12;
        nuctext.FontWeight = 'bold';
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

clearvars NuclearMeans NuclearSums;
disp('Collating results...')


    NuclearMeans(:,1) = Results_NuclearAnalysis.MeanCh1;
    NuclearMeans(:,2) = Results_NuclearAnalysis.MeanCh2;
    if Channels>2, NuclearMeans(:,3) = Results_NuclearAnalysis.MeanCh3; else end
    if Channels>3, NuclearMeans(:,4) = Results_NuclearAnalysis.MeanCh4; else end
    NuclearSums(:,1) = Results_NuclearAnalysis.SumCh1;
    NuclearSums(:,2) = Results_NuclearAnalysis.SumCh2;
    if Channels >2, NuclearSums(:,3) = Results_NuclearAnalysis.SumCh3; else end
    if Channels >3, NuclearSums(:,4) = Results_NuclearAnalysis.SumCh4; else end

if ZPlanes > 3 && (Results(f).ZFocus == 1 || Results(f).ZFocus == ZPlanes)
    Results(f).TotalNuclei = 0;
    Results(f).MeanIntensities = ""; %Each column is a channel.
    Results(f).SumIntensities = ""; %Each column is a channel.
else
    Results(f).TotalNuclei = sum(Results_NuclearAnalysis.NucleiNumber,'all');
    Results(f).MeanIntensities = NuclearMeans; %Each column is a channel.
    Results(f).SumIntensities = NuclearSums; %Each column is a channel.
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
[Results_Sorted,WellNumber] = WellSort(Channels,Results);
cd(Folder); cd Analysis;
save('AnalysisResultsSortedByWell.mat','Results_Sorted','-v7.3');

toc