%% Clear Workspace
clear all;
clc;

%% DIRECTORY STUFF (CHANGE THEM)

DAPI_inDir = 'Z:\Doug\Unsorted\09.13.18 E320A Stock Validation, Infectivity\5561\W1 - DAPI\*.tif'; %File path for DAPI images ('path\*.tif')
DAPI_Folder = 'Z:\Doug\Unsorted\09.13.18 E320A Stock Validation, Infectivity\5561\W1 - DAPI\'; % Folder path for DAPI images ('path\')

TAG_inDir = 'Z:\Doug\Unsorted\09.13.18 E320A Stock Validation, Infectivity\5561\W3 - VP1\*.tif'; %File path for TAg images ('path\*.tif')
TAG_Folder = 'Z:\Doug\Unsorted\09.13.18 E320A Stock Validation, Infectivity\5561\W3 - VP1\'; %Folder path for TAg images ('path\')

plate_rows = 8; %The number of rows of wells.
plate_columns = 12; %The number of columns of wells.
plate_locations = 5; %The number of images taken PER WELL
TAG_Threshold = 400; %THRESHOLD VALUE FOR MEAN TAG INTENSITY TO INDICATE "POSITIVE" INFECTED STATUS.
%% 

includeSubdirectories = true; %Tells the script to include subfolders within the input directories (inDir)

% Get relative path to all image files under inDir
DAPI_srcFiles = dir(DAPI_inDir); %Lists all files in the DAPI file directory indicated above
TAG_srcFiles = dir(TAG_inDir); %Lists all files in the TAg file directory indicated above
NTH = [1:plate_locations:(length(TAG_srcFiles))]; %Matrix for well summation and averaging (later)
Percent_infected = zeros(length(TAG_srcFiles),2); %Empty array that is populated later
TAG_meanaverage = zeros(1000,length(TAG_srcFiles));
figure('Position',[50 50 1200 800]) %[Left Bottom Width Height]
%% 
    for i = 1:length(DAPI_srcFiles) %Starts loop to analyze every DAPI image in folder using following commands
    DAPI_filename = strcat(DAPI_Folder,DAPI_srcFiles(i).name); %Opens image in file based on which loop the script is on
    TAG_filename = strcat(TAG_Folder,TAG_srcFiles(i).name);
    [pathstr,name,ext] = fileparts(TAG_filename); name02 = name(1:30);
    
    %%DAPI Nuclei Segmentation%%
    DAPI01 = imread(DAPI_filename); %Reads the image and names it "DAPI01."
    DAPI02 = im2double(DAPI01); %Converts DAPI01 to double format.
    
   LowBound = 20; UpBound = 250;
    
    DAPI02 = imadjust(DAPI02);
    bg = imopen(DAPI02,strel('disk',50));
    DAPI02 = DAPI02 - bg;
    BW = im2bw(DAPI02);
    BW02 = xor(bwareaopen(BW,LowBound), bwareaopen(BW,UpBound));
    BW03 = imerode(imdilate(imfill(BW02,'holes'),strel('disk',2)),strel('disk',2));
    DAPI_cc(i) = bwconncomp(BW03,8); %Identifies connected components (cells) within the binary image
    DAPI_cclabeled = labelmatrix(DAPI_cc(i));
    DAPI_ccRGB = label2rgb(DAPI_cclabeled,@autumn,'k','shuffle');
    
    %%TAG Filtering and Analysis%%
    TAG01 = imread(TAG_filename);
    TAG02 = im2double(TAG01);
    
    %Threshold and Size Filter%
   
    TAG_background = imopen(TAG01,strel('disk',100));
    TAG03 = TAG01 - TAG_background;
    TAG04 = TAG03>TAG_Threshold;
    TAG05 = TAG02(TAG04);
    %TAG_level01 = multithresh(TAG02,3);
    %TAG_bw01 = im2bw(TAG02,TAG_level01(1));
    %TAG_LowBound = 100; TAG_UpperBound = 700;
    %TAG_BW01 = xor(bwareaopen(TAG_bw01,TAG_LowBound), bwareaopen(TAG_bw01,TAG_UpperBound));
    %TAG_BW02 = imfill(TAG_BW01,'holes');
    %TAG03 = TAG02.*TAG_BW02;
    
    
    
    %%% Displays called images for each iteration (for monitoring)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clf
    
    TAG_PixelValues = regionprops('table',(DAPI_cc(i)),TAG03,'PixelValues');
    TAG_MeanIntensity = regionprops('table',(DAPI_cc(i)),TAG03,'mean');

    if DAPI_cc(i).NumObjects < 1
        TAG_meanarray = [0];
        TAG_meanaverage03 = TAG_meanarray;
    else
    
    
    TAG_meanarray = table2array(TAG_MeanIntensity); 
    TAG_meanaverage(1:length(TAG_meanarray),i) = TAG_meanarray(:,1);
    TAG_meanaverage02 = reshape(TAG_meanaverage,1,[]); TAG_meanaverage03 = TAG_meanaverage02(TAG_meanaverage02~=0);
    TAG_nucleiabovethresh = nnz(TAG_meanarray(:)>TAG_Threshold);
    DAPI_totalcells(i) = size(TAG_meanarray,1);
    Percent_infected(i,:) = [TAG_nucleiabovethresh (DAPI_totalcells(i))];
    
    %%For TAg Mean Intensity per cell%%
    TAG_infectedmeanfilter = TAG_meanaverage03>TAG_Threshold;
    TAG_uninfectedmeanfilter = TAG_meanaverage03<TAG_Threshold;
    TAG_infectedmeans = TAG_meanaverage03(TAG_infectedmeanfilter)';
    end
    
    Percent_monitoring = num2str(TAG_nucleiabovethresh/DAPI_cc(i).NumObjects*100);
    Cellcount_monitoring = num2str(DAPI_cc(i).NumObjects);
    
    progress_loop = num2str(i);
    progress_totalimages = num2str(plate_rows*plate_columns*plate_locations);
    progress_current = (i/(plate_rows*plate_columns*plate_locations))*10;
    progress_percent = num2str((i/(plate_rows*plate_columns*plate_locations))*100);
    nbins = 10;
    
    subplot(3,4,1);
        subimage(DAPI02.*3), title('DAPI'); axis off;
    subplot(3,4,2);
        subimage(DAPI_ccRGB), title('Segmented Nuclei'); axis off;
    subplot(3,4,[3,4]);
        h01 = histogram(DAPI_totalcells,nbins,'Normalization','probability');
        axis([0 500 0 0.5]);
        h01.EdgeColor = 'black';
        h01.BinWidth = 25;
        title('Cumulative Cell Number Distribution');
        text(25,0.4,Cellcount_monitoring,'Color','black','Fontsize',12);
        text(75,0.4,'cells in this image','Color','black','Fontsize',12);
    subplot(3,4,5);
        subimage(imadjust(TAG02)), title('TAg');  axis off;
    subplot(3,4,6);
        c01 = contourf(TAG_background(1:end,1:end),5);
        colormap hot; caxis([TAG_Threshold-.001 TAG_Threshold]); axis off;
        set(gca,'ydir','reverse'); title('TAg Background');
    subplot(3,4,7);
        subimage(imadjust(TAG03)), title('Backsubbed TAg'); axis off;
    subplot(3,4,8);
        subimage(imfuse(DAPI02.*3,TAG04)); title('Infected'); axis off;
    subplot(3,4,10); 
        h02 = histogram(TAG_meanarray,nbins,'Normalization','probability'); 
        h02.EdgeColor = 'black';
        h02.BinWidth = 100;
        axis([0 5000 0 0.5]); 
        title('TAg Intensity Distribution');
        line = rectangle('Position',[TAG_Threshold 0 50000-TAG_Threshold 1]);
        text(2500,0.4,Percent_monitoring,'Color','black','Fontsize',12);
        text(4100,0.4,'%','Color','black','Fontsize',12);
    subplot(3,4,[11,12]);
        h03 = histogram(TAG_meanaverage03,nbins,'Normalization','probability');
        axis([0 5000 0 0.5]);
        h03.EdgeColor = 'black';
        h03.BinWidth = 100;
        title('TAg Intensity Cumulative Distribution');
    subplot(3,4,9);
        file_monitoring = text(0,0,name02,'Color','black','Fontsize',10,'Interpreter','none');
        sub = rectangle('Position',[0 0 10 10],'LineStyle','none'); axis off; hold on;
        rec = rectangle('Position',[0 2 10 2],'LineStyle','-.'); axis off;
        prog = rectangle('Position',[0 2 progress_current 2],'FaceColor','red','LineStyle','-'); axis off; title('Progress');
        text(0,7,'% Complete:','Color','black','FontSize',12);
        text(6,7,progress_percent,'Color','black','FontSize',12);
        text(0,9,'Image','Color','black','FontSize',12);
        text(3,9,progress_loop,'Color','black','FontSize',12);
        text(5,9,'of','Color','black','FontSize',12);
        text(7,9,progress_totalimages,'Color','black','FontSize',12); drawnow;
    end
    
    WELL = zeros(96,1); %Another empty array that is populated later
    %% 
    
    for i = 1:(length(TAG_srcFiles)/plate_locations) %Starts a loop to determine the percent infected in each well
    WELL(i) = (sum(Percent_infected(NTH(:,i):(i*plate_locations),1)))/(sum(Percent_infected(NTH(:,i):(i*plate_locations),2)))*100; %Sums "count" and "CELLS" variables for each well and calculates percent infected
    end
    plate_length = plate_rows*plate_columns;
    TRIMMED = WELL(1:plate_length,1); %Removes zeros at end of WELL matrix.
    PLATE = reshape(TRIMMED,[plate_columns,plate_rows]); %Reshapes matrix into NxN matrix, as set by parameters at top of script. Yes, the COLUMNS and ROWS variables are reversed. It is because of the way the reshape command works.
    PLATE = PLATE'; %Flips the PLATE variable to create a matrix that should look just like your plate. You can copy/paste it directly into Excel.