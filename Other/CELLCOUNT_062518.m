%% Clear Workspace
clear all;
clc;

%% DIRECTORY STUFF (CHANGE THEM)

DAPI_inDir = 'C:\Users\keros\OneDrive\Documents\Plate images 7-5\5425 Plate One\DAPI 7-5\*.tif'; %File path for DAPI images ('path\*.tif')
DAPI_Folder = 'C:\Users\keros\OneDrive\Documents\Plate images 7-5\5425 Plate One\DAPI 7-5\'; % Folder path for DAPI images ('path\')

plate_locations = 3; %The number of images taken PER WELL
plate_rows = 8;
plate_columns = 12;

figure('Position',[100 100 1200 800])
%% 

includeSubdirectories = true; %Tells the script to include subfolders within the input directories (inDir)

% Get relative path to all image files under inDir
DAPI_srcFiles = dir(DAPI_inDir); %Lists all files in the DAPI file directory indicated above

%% 
    for i = 1:length(DAPI_srcFiles) %Starts loop to analyze every DAPI image in folder using following commands
    DAPI_filename = strcat(DAPI_Folder,DAPI_srcFiles(i).name); %Opens image in file based on which loop the script is on
    [pathstr,name,ext] = fileparts(DAPI_filename); name02 = name(1:30);
    
    DAPI01 = imread(DAPI_filename); %Calls the read image "DAPI"
    DAPI02 = im2double(DAPI01);
    %DAPI03 = imtophat(DAPI02,strel('disk',20));
    
    %Thresholding and Detection%
    %DAPI_level01 = multithresh(DAPI03,2); %Identifies the threshold value via Otsu's Method (2 layers)
    %DAPI_bw01 = im2bw(DAPI03,DAPI_level01(1)); %Converts DAPI image to binary using threshold determined in previous step
    %DAPI_bw02 = im2bw(DAPI03,DAPI_level01(2));
    %DAPI_BW01 = DAPI_bw01+DAPI_bw02;
    %DAPI_LowBound = 5; DAPI_UpperBound = 1100;
    %DAPI_BW02 = xor(bwareaopen(DAPI_BW01,DAPI_LowBound), bwareaopen(DAPI_BW01,DAPI_UpperBound)); %Removes outside of LowBound/UpperBound limits are removed
    %DAPI_BW03 = bwareaopen(imerode(DAPI_BW02,strel('disk',1)),5);
    %DAPI_cc(i) = bwconncomp(DAPI_BW03,8); %Identifies connected components (cells) within the binary image
    %DAPI_cclabeled = labelmatrix(DAPI_cc(i));
    %DAPI_ccRGB = label2rgb(DAPI_cclabeled,@autumn,'k','shuffle');
    %DAPI_BWoutline = imdilate(bwperim(DAPI_BW03),strel('disk',1));
    
    bg = imopen(DAPI02,strel('disk',100));
    DAPI03 = DAPI02 - bg;
    DAPI04 = imopen(DAPI03,strel('disk',5));
    BW = imbinarize(DAPI04);
    BW02 = bwareaopen(BW,300);
    BW03 = imdilate(imerode(BW02,strel('disk',3)),strel('disk',2));
    DAPI_cc(i) = bwconncomp(BW03,8); %Identifies connected components (cells) within the binary image
    DAPI_cclabeled = labelmatrix(DAPI_cc(i));
    DAPI_ccRGB = label2rgb(DAPI_cclabeled,@autumn,'k','shuffle');
    
    %Nuclear Area Calculation%
    Nuclei_area(i,:) = mean(table2array(regionprops('table',DAPI_cc(i),'area')));
    
    clf
    
    progress_loop = num2str(i);
    progress_totalimages = num2str(plate_rows*plate_columns*plate_locations);
    progress_current = (i/(plate_rows*plate_columns*plate_locations))*10;
    progress_percent = num2str((i/(plate_rows*plate_columns*plate_locations))*100);
    
    subplot(3,4,[1,2,5,6]); 
        subimage(imadjust(DAPI02)), title('Original'); axis off;
    subplot(3,4,[3,4,7,8]); 
        subimage(imfuse(imadjust(DAPI03),DAPI_BWoutline,'ColorChannels',[2 2 1])), title('Segmented Nuclei'); axis off;
    subplot(3,4,9);
        file_monitoring = text(0,0,name02,'Color','black','Fontsize',8,'Interpreter','none');
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
    
    for k = 1:length(DAPI_srcFiles)
        Nuclei_number(k) = DAPI_cc(k).NumObjects;
    end
   
    Nuclei_area_reshaped = reshape(Nuclei_area,[plate_locations,(plate_rows*plate_columns)]);
    Nuclei_number_reshaped = reshape(Nuclei_number,[plate_locations,(plate_rows*plate_columns)]);
    for m = 1:length(Nuclei_number_reshaped);
    Nuclei_number_wellsum(1,:) = sum(Nuclei_number_reshaped);
    Nuclei_area_wellaverage(1,:) = mean(Nuclei_area_reshaped);
    end
    Nuclei_number_sum_FINAL = reshape(Nuclei_number_wellsum,[plate_columns,plate_rows])';
    Nuclei_area_average_FINAL = reshape(Nuclei_area_wellaverage,[plate_columns,plate_rows])';
    
    beep on, beep, pause(1), beep;