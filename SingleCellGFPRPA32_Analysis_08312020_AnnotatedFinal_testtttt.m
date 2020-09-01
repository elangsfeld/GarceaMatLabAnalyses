asdasdgagad

%% README %%
%Creator: Douglas Peters, PhD
%Created: August 31st, 2020 (this version)
%Organization: Garcea Laboratory in the MCDB Department and 
%   BioFrontiers Institute, University of Colorado Boulder
%Purpose: This script was created to analyze changes in the localization of
%   fluorescent signal intensity. Specifically, this was created to measure
%   changes in GFP-RPA32 localization within cellular nuclei after
%   infection with murine polyomavirus (MuPyV). However, this script may be
%   useful for analysis of other fluorescent biomarkers in the nucleus.
%MatLab Version: Written in MATLAB R2020a.
%Prerequisites: Some MatLab add-ons may be required for this to run
%   properly. The image input was designed for .tif images, but may be
%   adaptable to other file types. No independent nuclear biomarkers is
%   required; this script identifies the nuclear border using the same
%   channel whose localization is being analyzed.
%
%Input: This script was designed to analyze single-cell vignettes (i.e., an
%   image series that features a single cell), saved as a .tif file. It may
%   be adaptable to track and analyze multiple cells, but that was beyond
%   the scope of this study.

%% Directory Information %%
clc; clear all; %Clears command window and deletes all existing variables in the workspace.
tic %Starts timer to report duration at the end of the script.

%%%User Input Required%%%
Folder = 'G:\DATA\11.25-28.16 Timecourse (Increasing NG59RA MOIs, mutants)\Vignettes\NG59 (5,11)\'; %Tells the script where to find the images that are to be analyzed. Must end with backslash.
inDir = append(Folder,'*.tif'); %Tells the script what file type to look for in the Folder path. This script has only been tested with tif files saved in FIJI.
ShowFigs = 0; %Choose to show figures during analysis (1=Yes,0=No). Turning this off makes the script run MUCH faster.
SaveFigs = 0; %Choose to save images of figures during analysis (1=Yes,0=No).
VRC_Multiplier = 1.5; %The mean nuclear intensity is multiplied by VRC_Multiplier to determine the floating threshold for identifying VRC pixels. Can be used to adjust VRC "stringency" (higher=more,lower=less).
%%%%%%%%%%%%%%%%%%%%%%%%%

srcFiles = dir(inDir); %Identifies the files within the inDir filepath.
if ShowFigs == 1, %Creates a figure to show the current image and analysis readings (only if ShowFigs ==1).
    figure('Position',[100 100 1000 600]) %Creates a figure with [Left Bottom Width Height] dimensions. Can be adjusted to fit screen size.
    nbins = 50; %Dictates the number of bins in the figure histogram. Can be adjusted to change granularity of histogram.
else end %Ends figure creation argument.

%% Nuclear Segmentation and Analysis %%

for v = 1:length(srcFiles), %Sets up a 'for' loop for each image series in Folder.
    Filename = strcat(Folder,srcFiles(v).name); %Reads the full filename.
    finfo = imfinfo(Filename); %Reads information about the file.
    NumFrames = size(finfo,1); %Reports the number of time frames in the image series.
    ResX = finfo(1).Width; %Reports the image resolution in the X dimension.
    ResY = finfo(1).Height; %Reports the image resolution in the Y dimension.
    
    clearvars RPAvals RPA_RP cc cc_VRC VRCareasum Readouts I I2 I3 Ibw Ibw2 Framecount I_VRC; %Clears these variables in between each image series to avoid overwrite issues.
    clf %Clears figure for next image series.

for f = 1:NumFrames %Sets up a 'for' loop for each frame in the current image series.
    I(:,:,f) = imread(Filename,f,'Info',finfo); %Reads the image based on the file info reported above.
    I2(:,:,f) = im2double(I(:,:,f)); %Converts the image to a double (0-1 scale).
    I3(:,:,f) = imopen(I2(:,:,f),strel('disk',(round(mean(ResX,ResY)/15)))); %Smoothes out high contrast areas in the image (e.g. high-intensity foci) to make segmentation more accurate.
    I3(:,:,f) = imlocalbrighten(I3(:,:,f)); %Brightens the diffuse nuclear regions of signal to improve segmentation.
    Ibw(:,:,f) = imbinarize(I3(:,:,f)); %Converts the image to binary (this is the segmentation step).
    Ibw2(:,:,f) = bwareaopen(Ibw(:,:,f),50); %Removes any groups of positive values (=1) in the binary image that are under 50 pixels in area.
    Ibw2(:,:,f) = imerode(imdilate(bwfill(Ibw2(:,:,f),'holes'),strel('disk',5)),strel('disk',5)); %These processing steps (fill, dilate, erode) smooth out the nuclear border.
    Ibw2(:,:,f) = bwareafilt(Ibw2(:,:,f),1); %Removes all but the largest group of +1 pixels from the image, which should be the nucleus.
    cc(f) = bwconncomp(Ibw2(:,:,f),8); %Makes the nuclear pixels a group that can be referenced later.
    Framecount(f,1) = f; %Creates an array for recording which frame is which.
    
    RPA_RP(f) = regionprops(cc(f),I2(:,:,f),'PixelValues','Circularity'); %Records the values of every nuclear pixel, as well as the circularity of the nuclear shape (1=perfect circle,0=line).
    RPAvals(f).PixelValues = RPA_RP(f).PixelValues; %Reports the pixels values to a different variable for easier access.
    RPAvals(f).NucCircularity = RPA_RP(f).Circularity; %Reports the nuclear circularity value to a different variable for easier access.
    RPAvals(f).NucMeanSignal = mean(RPAvals(f).PixelValues); %Calculates the mean values of all nuclear pixels.
    RPAvals(f).VRCSignalThresh = VRC_Multiplier*RPAvals(f).NucMeanSignal; %Calculates the floating threshold for identifying VRC pixels.
    RPAvals(f).VRCPixelidx = RPAvals(f).PixelValues>RPAvals(f).VRCSignalThresh; %Applies the floating VRC pixel threshold.
    RPAvals(f).VRCSignalsum = sum(RPAvals(f).PixelValues(RPAvals(f).VRCPixelidx)); %Calculates the integrated signal density of all VRC pixels.
    RPAvals(f).VRCSignalratio = mean((RPAvals(f).PixelValues(RPAvals(f).VRCPixelidx))./RPAvals(f).NucMeanSignal); %Calculates how much brighter VRC pixels are than the nuclear mean, on average across the nucleus.
    RPAvals(f).NucSignalsum = sum(RPAvals(f).PixelValues); %Calculates the integrated signal density of the whole nucleus.
    RPAvals(f).PrcntVRCSignal = RPAvals(f).VRCSignalsum/RPAvals(f).NucSignalsum*100; %Calculates the '% of Signal in VRCs'.
    
    I_VRC(:,:,f) = I2(:,:,f)>RPAvals(f).VRCSignalThresh; %Creates a binary image of all VRC pixels.
    cc_VRC(f) = bwconncomp(I_VRC(:,:,f),8); %Identifies connected +1 pixels into groups to get a VERY rough approximation of how many VRCs are present in the image.
    RPAvals(f).VRCnum = cc_VRC(f).NumObjects; %Reports the number of VRCs detected in the image. **NOT PRECISE**
    clearvars VRCarea NucPixVals; %Clears these variables between frames to avoid overwrite errors.
    if RPAvals(f).VRCnum > 0 %As long as there are VRC pixels in the image,
    for n = 1:RPAvals(f).VRCnum %Sets up a 'for' loop for each VRC detected.
        VRCareasum(n,1) = numel(cc_VRC(f).PixelIdxList{1,n}); %Reports the area (sq. pixels) of each VRC detected.
    end
    else VRCareasum(1,1) = 0; %Reports 0 as area value if there are no VRC pixels in the image.
    end
        
    RPAvals(f).VRCareasum = sum(VRCareasum); %Calculates the summed area of all VRC pixels.
    RPAvals(f).PrcntVRCarea = RPAvals(f).VRCareasum/numel(RPAvals(f).PixelValues)*100; %Calculates the '% of Nuclear Area in VRCs'.
    RPAvals(f).VRCareamean = mean(VRCareasum); %Calculates the mean area of VRCs detected.  **NOT PRECISE**
 
    Readouts(f,1) = RPAvals(f).NucCircularity; %Moves circularity value to a new variable.
    Readouts(f,2) = RPAvals(f).NucMeanSignal; %Moves nuclear mean intensity value to a new variable.
    Readouts(f,3) = RPAvals(f).VRCSignalratio; %Moves VRC:Nuclear signal ratio value to a new variable.
    Readouts(f,4) = RPAvals(f).PrcntVRCSignal; %Moves % of Signal in VRCs value to a new variable.
    Readouts(f,5) = RPAvals(f).VRCnum; %Moves VRC number value to a new variable.  **NOT PRECISE**
    Readouts(f,6) = RPAvals(f).PrcntVRCarea; %Moves % of Nuclear Area in VRCS value to a new variable.
    Readouts(f,7) = RPAvals(f).VRCareamean; %Moves mena VRC area value to a new variable.  **NOT PRECISE**
    
if ShowFigs == 1   %If you want the script to print the figures to the screen....
    subplot(3,5,[1,2,6,7]), subimage(imfuse(imadjust(I(:,:,f)),bwperim(Ibw2(:,:,f),8),'ColorChannels',[2 1 2])), title('Nuclear Segmentation'); axis off; hold on;
        text(2,3,srcFiles(v).name,'Color',[1 1 1]); hold off; %Reports the nuclear segmentation of the current image.
    subplot(3,5,[3,4,5]), his = histogram(RPAvals(f).PixelValues,nbins,'Normalization','probability');
        axis([0 1 0 0.25]), title('Pixel Int. Distr. w/in Nuc'), his.BinWidth = 0.01; %Reports the signal intensity histogram of the current image.
    subplot(3,5,8), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).PrcntVRCSignal,'filled'); grid on; title('%Signal w/in VRCs'); xlabel('Frame'); hold off; %Reports the % of Signal in VRCs for each frame analyzed so far.
    subplot(3,5,9), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).PrcntVRCarea,'filled'); grid on; title('%Area w/in VRCs'); xlabel('Frame'); hold off; %Reports the % of Nuclear area in VRCs for each frame analyzed so far.
    subplot(3,5,10), axis([1 NumFrames 0 1.1]); hold on;
        scatter(f,RPAvals(f).NucCircularity,'filled'); grid on; title('Nuc. Circularity'); xlabel('Frame'); hold off; %Reports the Nuclear Circularity value for each frame analyzed so far.
    subplot(3,5,11), axis([1 NumFrames 0 5000]); hold on;
        scatter(f,numel(RPAvals(f).PixelValues),'filled'); grid on; title('Total Nuclear Area'); xlabel('Frame'); hold off; %Reports the Nuclear Area value for each frame analyzed so far.
    subplot(3,5,12), axis([1 NumFrames 0 4]); hold on;
        scatter(f,RPAvals(f).VRCSignalratio,'filled'); grid on; title('VRC:Nuc Ratio'); xlabel('Frame'); hold off; %Reports the VRC:Nuc Signal ratio for each frame analyzed so far.
    subplot(3,5,13), axis([1 NumFrames 0 0.5]); hold on;
        scatter(f,RPAvals(f).NucMeanSignal,'filled'); grid on; title('Mean Nuclear Signal'); xlabel('Frame'); hold off; %Reports the mean nuclear signal intensity value for each frame analyzed so far.
    subplot(3,5,14), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).VRCnum,'filled'); grid on; title('VRC number'); xlabel('Frame'); hold off; %Reports the number of VRCs detected.  **NOT PRECISE**
    subplot(3,5,15), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).VRCareamean,'filled'); grid on; title('Mean VRC Area'); xlabel('Frame'); hold off; drawnow; %Reports the mean VRC area of those detected.  **NOT PRECISE**
else
end
if f == NumFrames && ShowFigs ==1 %If this is the last image (i.e., frame) of the series and you chose to show the analysis figure...
        subplot(3,4,[3,4]), his = histogram(double(RPAvals(f).PixelValues),nbins,'Normalization','probability'); hold on; %This will plot the histogram of the first frame against that of the final frame.
        his_og = histogram(RPAvals(1).PixelValues,nbins,'Normalization','probability');
        axis([0 1 0 0.25]); title('Pixel Int. Distr. w/in Nuc');
        his.BinWidth = 0.01; his_og.BinWidth = 0.01;
        legend ('Last','First'); drawnow;
end

if f == NumFrames %If this is the last image (i.e., frame) of the series...
    cd(Folder); %Changes directory to Folder.
    figfilename = strcat(srcFiles(v).name(1:end-4),' Analysis.pdf'); %Creates the filename that will be used to record an image of the analysis figure (as a pdf).
    matfilename = strcat(srcFiles(v).name(1:end-4),' Workspace.mat'); %Creates the filename that will be used to record the entire workspace of variables.
if v == 1, %If this is the first image series being analyzed...
    mkdir Analysis; %Makes a folder in Folder named 'Analysis'.
    cd('Analysis'); %Changes directory to newly-created Analysis folder.
else cd('Analysis'); %If this isn't the first image being analyzed, then there should already be an Analysis folder, which is moved to by this command.
end 
save(matfilename); %Saves the workspace of this cell's analysis, using the filename created above.
if ShowFigs == 1 && SaveFigs == 1 %If you want to show and save the figure images...
    exportgraphics(gcf,figfilename); %Saves the figure image as a pdf using the filename created above.
else
end
end

end %Ends the analysis for this image series.

%%% The following lines create the Results variable, which compiles these
%%% analysis values for all the image series that were analyzed.
Results(v).NucCircularity = num2cell(Readouts(:,1)); %Nuclear circularity
Results(v).NucMeanSignal = num2cell(Readouts(:,2)); %Mean nuclear signal intensity
Results(v).VRCSignalratio = num2cell(Readouts(:,3)); %Mean VRC:Nuc signal ratio
Results(v).PrcntVRCSignal = num2cell(Readouts(:,4)); % '% of Signal in VRCs'
Results(v).VRCnum = num2cell(Readouts(:,5)); %Number of VRCs detected **NOT PRECISE**
Results(v).PrcntVRCarea = num2cell(Readouts(:,6)); % '% of Nuclear Area in VRCs'
Results(v).VRCareamean = num2cell(Readouts(:,7)); % 'Mean area of VRCs detected **NOT PRECISE**
end

save('Full Results.mat','Results'); %Saves the Results variable to a .mat file in the Analysis folder created above.
close all; %Closes all figures that are open.
toc %Ends timer for script and prints result on command line.