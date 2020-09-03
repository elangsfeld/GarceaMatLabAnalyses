blahblahb

%% Directory Info %%
clc; clear all; %Clears command window and deletes all existing variables in the workspace.
tic

inDir = 'G:\DATA\11.25-28.16 Timecourse (Increasing NG59RA MOIs, mutants)\Vignettes\NG59 (5,11)\*.tif';
Folder = 'G:\DATA\11.25-28.16 Timecourse (Increasing NG59RA MOIs, mutants)\Vignettes\NG59 (5,11)\';
srcFiles = dir(inDir);

ShowandSaveFigs = 1;

VRC_Multiplier = 1.5; %The mean nuclear intensity is multiplied by VRC_Multiplier to determine the floating threshold for identifying VRC pixels.

if ShowandSaveFigs == 1, 
    figure('Position',[100 100 1000 600]) %[Left Bottom Width Height]
    nbins = 50; %Dictates the number of bins in the figure histogram.
else end

%% Nuclei Segmentation %%

for v = 1:length(srcFiles),
    Filename = strcat(Folder,srcFiles(v).name);
    finfo = imfinfo(Filename);
    NumFrames = size(finfo,1);
    ResX = finfo(1).Width;
    ResY = finfo(1).Height;
    
    clearvars RPAvals RPA_RP cc cc_VRC VRCareasum Readouts I I2 I3 Ibw Ibw2 Framecount I_VRC;
    clf

for f = 1:NumFrames
    I(:,:,f) = imread(Filename,f,'Info',finfo);
    I2(:,:,f) = im2double(I(:,:,f));
    I3(:,:,f) = imopen(I2(:,:,f),strel('disk',(round(mean(ResX,ResY)/15))));
    I3(:,:,f) = imlocalbrighten(I3(:,:,f));
    Ibw(:,:,f) = imbinarize(I3(:,:,f));
    Ibw2(:,:,f) = bwareaopen(Ibw(:,:,f),50);
    Ibw2(:,:,f) = imerode(imdilate(bwfill(Ibw2(:,:,f),'holes'),strel('disk',5)),strel('disk',5));
    Ibw2(:,:,f) = bwareafilt(Ibw2(:,:,f),1);
    cc(f) = bwconncomp(Ibw2(:,:,f),8);
    Framecount(f,1) = f;
    
    RPA_RP(f) = regionprops(cc(f),I2(:,:,f),'PixelValues','Circularity');
    RPAvals(f).PixelValues = RPA_RP(f).PixelValues;
    RPAvals(f).NucCircularity = RPA_RP(f).Circularity;
    RPAvals(f).NucMeanSignal = mean(RPAvals(f).PixelValues);
    RPAvals(f).VRCSignalThresh = VRC_Multiplier*RPAvals(f).NucMeanSignal;
    RPAvals(f).VRCPixelidx = RPAvals(f).PixelValues>RPAvals(f).VRCSignalThresh;
    RPAvals(f).VRCSignalsum = sum(RPAvals(f).PixelValues(RPAvals(f).VRCPixelidx));
    RPAvals(f).VRCSignalratio = mean((RPAvals(f).PixelValues(RPAvals(f).VRCPixelidx))./RPAvals(f).NucMeanSignal);
    RPAvals(f).NucSignalsum = sum(RPAvals(f).PixelValues);
    RPAvals(f).PrcntVRCSignal = RPAvals(f).VRCSignalsum/RPAvals(f).NucSignalsum*100;
    
    
    
    I_VRC(:,:,f) = I2(:,:,f)>RPAvals(f).VRCSignalThresh;
    cc_VRC(f) = bwconncomp(I_VRC(:,:,f),8);
    
    RPAvals(f).VRCnum = cc_VRC(f).NumObjects;
    
    clearvars VRCarea NucPixVals;
    if RPAvals(f).VRCnum > 0
    for n = 1:RPAvals(f).VRCnum
        VRCareasum(n,1) = numel(cc_VRC(f).PixelIdxList{1,n});
    end
    else VRCareasum(1,1) = 0;
    end
        
    RPAvals(f).VRCareasum = sum(VRCareasum);
    RPAvals(f).PrcntVRCarea = RPAvals(f).VRCareasum/numel(RPAvals(f).PixelValues)*100;
    RPAvals(f).VRCareamean = mean(VRCareasum);
 
    Readouts(f,1) = RPAvals(f).NucCircularity;
    Readouts(f,2) = RPAvals(f).NucMeanSignal;
    Readouts(f,3) = RPAvals(f).VRCSignalratio;
    Readouts(f,4) = RPAvals(f).PrcntVRCSignal;
    Readouts(f,5) = RPAvals(f).VRCnum;
    Readouts(f,6) = RPAvals(f).PrcntVRCarea;
    Readouts(f,7) = RPAvals(f).VRCareamean;    
    
if ShowandSaveFigs == 1   
    %clf %Clears figure from the previous image.
    subplot(3,5,[1,2,6,7]), subimage(imfuse(imadjust(I(:,:,f)),bwperim(Ibw2(:,:,f),8),'ColorChannels',[2 1 2])), title('Nuclear Segmentation'); axis off; hold on;
        text(2,3,srcFiles(v).name,'Color',[1 1 1]); hold off;
    subplot(3,5,[3,4,5]), his = histogram(RPAvals(f).PixelValues,nbins,'Normalization','probability');
        axis([0 1 0 0.25]), title('Pixel Int. Distr. w/in Nuc'), his.BinWidth = 0.01;
    subplot(3,5,8), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).PrcntVRCSignal,'filled'); grid on; title('%Signal w/in VRCs'); xlabel('Frame'); hold off; 
    subplot(3,5,9), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).PrcntVRCarea,'filled'); grid on; title('%Area w/in VRCs'); xlabel('Frame'); hold off; 
    subplot(3,5,10), axis([1 NumFrames 0 1.1]); hold on;
        scatter(f,RPAvals(f).NucCircularity,'filled'); grid on; title('Nuc. Circularity'); xlabel('Frame'); hold off; 
    subplot(3,5,11), axis([1 NumFrames 0 5000]); hold on;
        scatter(f,numel(RPAvals(f).PixelValues),'filled'); grid on; title('Total Nuclear Area'); xlabel('Frame'); hold off; 
    subplot(3,5,12), axis([1 NumFrames 0 4]); hold on;
        scatter(f,RPAvals(f).VRCSignalratio,'filled'); grid on; title('VRC:Nuc Ratio'); xlabel('Frame'); hold off; 
    subplot(3,5,13), axis([1 NumFrames 0 0.5]); hold on;
        scatter(f,RPAvals(f).NucMeanSignal,'filled'); grid on; title('Mean Nuclear Signal'); xlabel('Frame'); hold off;
    subplot(3,5,14), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).VRCnum,'filled'); grid on; title('VRC number'); xlabel('Frame'); hold off;
    subplot(3,5,15), axis([1 NumFrames 0 50]); hold on;
        scatter(f,RPAvals(f).VRCareamean,'filled'); grid on; title('Mean VRC Area'); xlabel('Frame'); hold off; drawnow;
    

    if f == NumFrames,
        
        subplot(3,4,[3,4]), his = histogram(double(RPAvals(f).PixelValues),nbins,'Normalization','probability'); hold on;
        his_og = histogram(RPAvals(1).PixelValues,nbins,'Normalization','probability');
        axis([0 1 0 0.25]); title('Pixel Int. Distr. w/in Nuc');
        his.BinWidth = 0.01; his_og.BinWidth = 0.01;
        legend ('Last','First'); drawnow;
        
    cd(Folder);
    figfilename = strcat(srcFiles(v).name(1:end-4),' Analysis.pdf');
    matfilename = strcat(srcFiles(v).name(1:end-4),' Workspace.mat');
        
    if v == 1,
    mkdir Analysis;
    cd('Analysis'); 
    else cd('Analysis'); end
    exportgraphics(gcf,figfilename);
    save(matfilename);
   
    else
    end
else
end
end

Results(v).NucCircularity = num2cell(Readouts(:,1));
Results(v).NucMeanSignal = num2cell(Readouts(:,2));
Results(v).VRCSignalratio = num2cell(Readouts(:,3));
Results(v).PrcntVRCSignal = num2cell(Readouts(:,4));
Results(v).VRCnum = num2cell(Readouts(:,5));
Results(v).PrcntVRCarea = num2cell(Readouts(:,6));
Results(v).VRCareamean = num2cell(Readouts(:,7));
    %pause(0.5)
end

save('Full Results.mat','Results');

toc
