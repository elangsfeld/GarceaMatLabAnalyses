clear all; clc;
tic
cd(userpath);

inDir = 'C:\Users\dplad\Desktop\SDC test\40x\*.nd2';
Folder = 'C:\Users\dplad\Desktop\SDC test\40x\';

srcFiles = dir(inDir);

%%% Requires INPUT %%%
LPercentileThresh(1,1:4) = [99.5 99.75 99.75 99.75];
UPercentileThresh(1,1:4) = [100 100 100 100];

Channels = 2;
drawVRCimages = 1; %Option to draw and save VRC boxes on image as PDF (1=yes, 0=no).

DBSCAN_Rad(1,1:4) = [2 2 2 2];
DBSCAN_MinPts(1,1:4) = [5 5 5 5];

RS_Rad(1,1:4) = [5 5 5 5];

VRC_ID = [1 0 0 0]; %Chooses which channel(s) to use for VRC identification (can be up to 2).
VRC_IDCH = 1; %Chooses which channel to use for VRC point seeding (for filtering rangesearch).

AnalysisVolume = 0;
AnalysisPCC = 0;
%%% %%%

VRC_ID_low = 'Must choose at least one channel for VRC identification (VRC_ID).';
VRC_ID_high = 'Must choose UP TO 2 channels for VRC identification (VRC_ID).';
varNames = {'MinY','MaxY','MinX','MaxX','CenY','CenX'};
varTypes = {'double','double','double','double','double','double'};
varNames_RS = {'Ch1MinY','Ch1MaxY','Ch1MinX','Ch1MaxX','Ch2MinY','Ch2MaxY','Ch2MinX','Ch2MaxX','Ch3MinY','Ch3MaxY','Ch3MinX','Ch3MaxX','Ch4MinY','Ch4MaxY','Ch4MinX','Ch4MaxX'};
varTypes_RS = {'double','double','double','double','double','double','double','double','double','double','double','double','double','double','double','double'};
varNames_VRC = {'VRC_MinY','VRC_MaxY','VRC_MinX','VRC_MaxX'};
varTypes_VRC = {'double','double','double','double'};
varNames_Volumes = {'Ch1Pix','Ch2Pix','Ch3Pix','Ch4Pix'};
varTypes_Volumes = {'double','double','double','double'};
varNames_PCCs = {'1-2','1-3','1-4','2-3','2-4','3-4'};
varTypes_PCCs = {'double','double','double','double','double','double'};
VRC_Idx = [1,5,9,13; 2,6,10,14; 3,7,11,15; 4,8,12,16];
VRC_ID2 = [VRC_ID;VRC_ID;VRC_ID;VRC_ID];
VRC_Idx2 = VRC_Idx(logical(VRC_ID2)); VRC_Idx2 = reshape(VRC_Idx2,[4 sum(VRC_ID(:))]);
if sum(VRC_ID(:)) < 2,
    VRC_Idx2(1:4,2) = [NaN;NaN;NaN;NaN];
else end


for f = 1:length(srcFiles)
    progress = (f/length(srcFiles))*100;
    progress2 = sprintf('Analyzing image %d of %d; %0.2f%c complete.',f,length(srcFiles),progress,'%');
    disp(progress2)
    
    filename = strcat(Folder,srcFiles(f).name);
    
    cd(Folder);
    
    I = bfopen(filename);
    
    ResX = size(I{1,1}{1,1},1);
    ResY = size(I{1,1}{1,1},2);
    Slices = (length(I{1,1})/Channels);
    Blank3D = zeros(ResX,ResY,Slices);
    Blank2D = zeros(ResX,ResY);
    Planes = Channels*Slices;
    
    VRCs = {}; %Creates cell variable for population later. Columns refer to each fluorescent channel.
    
    Ch1 = uint16(Blank3D);
    Ch2 = uint16(Blank3D);
    if Channels>2, Ch3 = uint16(Blank3D); else end;
    if Channels>3, Ch4 = uint16(Blank3D); else end;

%% SEPARATES EACH FLUORESCENT CHANNEL INTO ITS OWN STACK OF IMAGES %%
for i = 1:Slices
    Ch1_planes(i,1) = 1+(Channels*i-Channels);
    Ch2_planes(i,1) = 2+(Channels*i-Channels);
    if Channels>2, Ch3_planes(i,1) = 3+(Channels*i-Channels); else end
    if Channels>3, Ch4_planes(i,1) = 4+(Channels*i-Channels); else end
end

for m = 1:Slices
    Ch1(:,:,m) = I{1,1}{Ch1_planes(m,1),1};
    Ch1_IntDen_Slice(f,m) = sum(sum(Ch1(:,:,m)));
    [Ch1_Max(f,1), Ch1_MaxSlice(f,1)] = max(Ch1_IntDen_Slice(f,:));
    
    Ch2(:,:,m) = I{1,1}{Ch2_planes(m,1),1};
    Ch2_IntDen_Slice(f,m) = sum(sum(Ch2(:,:,m)));
    [Ch2_Max(f,1), Ch2_MaxSlice(f,1)] = max(Ch2_IntDen_Slice(f,:));
    
    if Channels>2, Ch3(:,:,m) = I{1,1}{Ch3_planes(m,1),1};
    Ch3_IntDen_Slice(f,m) = sum(sum(Ch3(:,:,m)));
    [Ch3_Max(f,1), Ch3_MaxSlice(f,1)] = max(Ch3_IntDen_Slice(f,:));else end
    
    if Channels>3, Ch4(:,:,m) = I{1,1}{Ch4_planes(m,1),1}; 
    Ch4_IntDen_Slice(f,m) = sum(sum(Ch4(:,:,m)));
    [Ch4_Max(f,1), Ch4_MaxSlice(f,1)] = max(Ch4_IntDen_Slice(f,:));else end   
end
%% THRESHOLDING BASED ON PERCENTILES ENTERED AT TOP OF SCRIPT %%
LCh1_thresh = prctile(Ch1,LPercentileThresh(1,1),'all');
UCh1_thresh = prctile(Ch1,UPercentileThresh(1,1),'all');
Ch1_BWMask = (UCh1_thresh>=Ch1) & (Ch1>LCh1_thresh); 
Ch1_seg = uint16(double(Ch1).*Ch1_BWMask);

LCh2_thresh = prctile(Ch2,LPercentileThresh(1,2),'all');
UCh2_thresh = prctile(Ch2,UPercentileThresh(1,2),'all');
Ch2_BWMask = (UCh2_thresh>=Ch2) & (Ch2>LCh2_thresh); 
Ch2_seg = uint16(double(Ch2).*Ch2_BWMask);

if Channels>2
    LCh3_thresh = prctile(Ch3,LPercentileThresh(1,3),'all');
    UCh3_thresh = prctile(Ch3,UPercentileThresh(1,3),'all');
    Ch3_BWMask = (UCh3_thresh>=Ch3) & (Ch3>LCh3_thresh); 
    Ch3_seg = uint16(double(Ch3).*Ch3_BWMask);
else end

if Channels>3
    LCh4_thresh = prctile(Ch4,LPercentileThresh(1,4),'all');
    UCh4_thresh = prctile(Ch4,UPercentileThresh(1,4),'all');
    Ch4_BWMask = (UCh4_thresh>=Ch4) & (Ch4>LCh4_thresh); 
    Ch4_seg = uint16(double(Ch4).*Ch4_BWMask);
else end
%% XY COORDINATE LISTS OF POSITIVE PIXELS FROM BINARY MASKS, THEN DBSCAN ANALYSIS%%

if VRC_ID(1,1) == 1,
clearvars Ch1clusterROW Ch1clusterCOL Ch1clusterXY Ch1_DBSCAN Ch1_DBSCANfilter Ch1_DBSCAN2 Ch1_DBSCAN2_Max Ch1clusterROW2 Ch1clusterCOL2 Ch1clusterXY2 Ch1_Boxes Ch1_Cens; 
[Ch1clusterROW,Ch1clusterCOL] = find(Ch1_BWMask(:,:,Ch1_MaxSlice(f,1))); %IDs above-threshold pixels in the brightest z-slice.
Ch1clusterXY = [Ch1clusterCOL, Ch1clusterROW]; %Combines positive coordinates in one array.
Ch1_DBSCAN = dbscan(Ch1clusterXY,DBSCAN_Rad(1,1),DBSCAN_MinPts(1,1)); %Performs DBSCAN analysis based on criteria in the top of script.
Ch1_DBSCANfilter = Ch1_DBSCAN>0; %Makes logical for DBSCAN groups (>0 ID). Outliers get -1 ID that needs to be filtered out.
Ch1_DBSCAN2 = Ch1_DBSCAN(Ch1_DBSCANfilter); %Filters OG DBSCAN to remove -1 IDs.
Ch1_DBSCAN2_Max = max(Ch1_DBSCAN2(:)); %Finds max value in DBSCAN to determine how many clusters were detected.
Ch1clusterROW2 = Ch1clusterROW(Ch1_DBSCANfilter); %Filters rows
Ch1clusterCOL2 = Ch1clusterCOL(Ch1_DBSCANfilter); %and columns
Ch1clusterXY2 = [Ch1clusterCOL2, Ch1clusterROW2]; %and puts them into array.

Ch1_Boxes = table('Size',[Ch1_DBSCAN2_Max length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames); %Creates table variable for bounding box information.
Ch1_gscatter = gscatter(Ch1clusterXY2(:,1),Ch1clusterXY2(:,2),Ch1_DBSCAN2,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.

for b = 1:Ch1_DBSCAN2_Max %This loops finds the min/max Y/X values for each "group IDed by gscatter function above, then calculates their centroids.
Ch1_Boxes.MinY(b) = min(Ch1_gscatter(b,1).YData);
Ch1_Boxes.MaxY(b) = max(Ch1_gscatter(b,1).YData);
Ch1_Boxes.MinX(b) = min(Ch1_gscatter(b,1).XData);
Ch1_Boxes.MaxX(b) = max(Ch1_gscatter(b,1).XData);

Ch1_Boxes.CenY(b) = mean(Ch1_gscatter(b,1).YData);
Ch1_Boxes.CenX(b) = mean(Ch1_gscatter(b,1).XData);
end

Ch1_Cens(1:length(Ch1_Boxes.CenY),1) = Ch1_Boxes.CenY(:);
Ch1_Cens(1:length(Ch1_Boxes.CenX),2) = Ch1_Boxes.CenX(:);


else end

if VRC_ID(1,2) == 1,
    clearvars Ch2clusterROW Ch2clusterCOL Ch2clusterXY Ch2_DBSCAN Ch2_DBSCANfilter Ch2_DBSCAN2 Ch2_DBSCAN2_Max Ch2clusterROW2 Ch2clusterCOL2 Ch2clusterXY2  Ch2_Boxes Ch2_Cens;
    [Ch2clusterROW,Ch2clusterCOL] = find(Ch2_BWMask(:,:,Ch2_MaxSlice(f,1))); %IDs above-threshold pixels in the brightest z-slice.
    Ch2clusterXY = [Ch2clusterCOL, Ch2clusterROW]; %Combines positive coordinates in one array.
    Ch2_DBSCAN = dbscan(Ch2clusterXY,DBSCAN_Rad(1,2),DBSCAN_MinPts(1,2)); %Performs DBSCAN analysis based on criteria in the top of script.
    Ch2_DBSCANfilter = Ch2_DBSCAN>0; %Makes logical for DBSCAN groups (>0 ID). Outliers get -1 ID that needs to be filtered out.
    Ch2_DBSCAN2 = Ch2_DBSCAN(Ch2_DBSCANfilter); %Filters OG DBSCAN to remove -1 IDs.
    Ch2_DBSCAN2_Max = max(Ch2_DBSCAN2(:)); %Finds max value in DBSCAN to determine how many clusters were detected.
    Ch2clusterROW2 = Ch2clusterROW(Ch2_DBSCANfilter); %Filters rows
    Ch2clusterCOL2 = Ch2clusterCOL(Ch2_DBSCANfilter); %and columns
    Ch2clusterXY2 = [Ch2clusterCOL2, Ch2clusterROW2]; %and puts them into array.

    Ch2_Boxes = table('Size',[Ch2_DBSCAN2_Max length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames); %Creates table variable for bounding box information.
    Ch2_gscatter = gscatter(Ch2clusterXY2(:,1),Ch2clusterXY2(:,2),Ch2_DBSCAN2,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.

    for b = 1:Ch2_DBSCAN2_Max %This loops finds the min/max Y/X values for each "group IDed by gscatter function above, then calculates their centroids.
    Ch2_Boxes.MinY(b) = min(Ch2_gscatter(b,1).YData);
    Ch2_Boxes.MaxY(b) = max(Ch2_gscatter(b,1).YData);
    Ch2_Boxes.MinX(b) = min(Ch2_gscatter(b,1).XData);
    Ch2_Boxes.MaxX(b) = max(Ch2_gscatter(b,1).XData);
    
    Ch2_Boxes.CenY(b) = mean(Ch2_gscatter(b,1).YData);
    Ch2_Boxes.CenX(b) = mean(Ch2_gscatter(b,1).XData);
    end
        
    Ch2_Cens(1:length(Ch2_Boxes.CenY),1) = Ch2_Boxes.CenY(:);
    Ch2_Cens(1:length(Ch2_Boxes.CenX),2) = Ch2_Boxes.CenX(:);
    
else end

if Channels>2 && VRC_ID(1,3) == 1,
    clearvars Ch3clusterROW Ch3clusterCOL Ch3clusterXY Ch3_DBSCAN Ch3_DBSCANfilter Ch3_DBSCAN2 Ch3_DBSCAN2_Max Ch3clusterROW2 Ch3clusterCOL2 Ch3clusterXY2 Ch3_Boxes Ch3_Cens;
    [Ch3clusterROW,Ch3clusterCOL] = find(Ch3_BWMask(:,:,Ch3_MaxSlice(f,1))); %IDs above-threshold pixels in the brightest z-slice.
    Ch3clusterXY = [Ch3clusterCOL, Ch3clusterROW]; %Combines positive coordinates in one array.
    Ch3_DBSCAN = dbscan(Ch3clusterXY,DBSCAN_Rad(1,3),DBSCAN_MinPts(1,3)); %Performs DBSCAN analysis based on criteria in the top of script.
    Ch3_DBSCANfilter = Ch3_DBSCAN>0; %Makes logical for DBSCAN groups (>0 ID). Outliers get -1 ID that needs to be filtered out.
    Ch3_DBSCAN2 = Ch3_DBSCAN(Ch3_DBSCANfilter); %Filters OG DBSCAN to remove -1 IDs.
    Ch3_DBSCAN2_Max = max(Ch3_DBSCAN2(:)); %Finds max value in DBSCAN to determine how many clusters were detected.
    Ch3clusterROW2 = Ch3clusterROW(Ch3_DBSCANfilter); %Filters rows
    Ch3clusterCOL2 = Ch3clusterCOL(Ch3_DBSCANfilter); %and columns
    Ch3clusterXY2 = [Ch3clusterCOL2, Ch3clusterROW2]; %and puts them into array.
    
    Ch3_Boxes = table('Size',[Ch3_DBSCAN2_Max length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames); %Creates table variable for bounding box information.
    Ch3_gscatter = gscatter(Ch3clusterXY2(:,1),Ch3clusterXY2(:,2),Ch3_DBSCAN2,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.

    for b = 1:Ch3_DBSCAN2_Max %This loops finds the min/max Y/X values for each "group IDed by gscatter function above, then calculates their centroids.
    Ch3_Boxes.MinY(b) = min(Ch3_gscatter(b,1).YData);
    Ch3_Boxes.MaxY(b) = max(Ch3_gscatter(b,1).YData);
    Ch3_Boxes.MinX(b) = min(Ch3_gscatter(b,1).XData);
    Ch3_Boxes.MaxX(b) = max(Ch3_gscatter(b,1).XData);
    
    Ch3_Boxes.CenY(b) = mean(Ch3_gscatter(b,1).YData);
    Ch3_Boxes.CenX(b) = mean(Ch3_gscatter(b,1).XData);
    end
    
    Ch3_Cens(1:length(Ch3_Boxes.CenY),1) = Ch3_Boxes.CenY(:);
    Ch3_Cens(1:length(Ch3_Boxes.CenX),2) = Ch3_Boxes.CenX(:);

else end

if Channels>3 && VRC_ID(1,4) == 1,
    clearvars Ch4clusterROW Ch4clusterCOL Ch4clusterXY Ch4_DBSCAN Ch4_DBSCANfilter Ch4_DBSCAN2 Ch4_DBSCAN2_Max Ch4clusterROW2 Ch4clusterCOL2 Ch4clusterXY2 Ch4_Boxes Ch4_Cens;
    [Ch4clusterROW,Ch4clusterCOL] = find(Ch4_BWMask(:,:,Ch4_MaxSlice(f,1))); %IDs above-threshold pixels in the brightest z-slice.
    Ch4clusterXY = [Ch4clusterCOL, Ch4clusterROW]; %Combines positive coordinates in one array.
    Ch4_DBSCAN = dbscan(Ch4clusterXY,DBSCAN_Rad(1,4),DBSCAN_MinPts(1,4)); %Performs DBSCAN analysis based on criteria in the top of script.
    Ch4_DBSCANfilter = Ch4_DBSCAN>0; %Makes logical for DBSCAN groups (>0 ID). Outliers get -1 ID that needs to be filtered out.
    Ch4_DBSCAN2 = Ch4_DBSCAN(Ch4_DBSCANfilter); %Filters OG DBSCAN to remove -1 IDs.
    Ch4_DBSCAN2_Max = max(Ch4_DBSCAN2(:)); %Finds max value in DBSCAN to determine how many clusters were detected.
    Ch4clusterROW2 = Ch4clusterROW(Ch4_DBSCANfilter); %Filters rows
    Ch4clusterCOL2 = Ch4clusterCOL(Ch4_DBSCANfilter); %and columns
    Ch4clusterXY2 = [Ch4clusterCOL2, Ch4clusterROW2]; %and puts them into array.
    
    Ch4_Boxes = table('Size',[Ch4_DBSCAN2_Max length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames); %Creates table variable for bounding box information.
    Ch4_gscatter = gscatter(Ch4clusterXY2(:,1),Ch4clusterXY2(:,2),Ch4_DBSCAN2,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.

    for b = 1:Ch4_DBSCAN2_Max %This loops finds the min/max Y/X values for each "group IDed by gscatter function above, then calculates their centroids.
    Ch4_Boxes.MinY(b) = min(Ch4_gscatter(b,1).YData);
    Ch4_Boxes.MaxY(b) = max(Ch4_gscatter(b,1).YData);
    Ch4_Boxes.MinX(b) = min(Ch4_gscatter(b,1).XData);
    Ch4_Boxes.MaxX(b) = max(Ch4_gscatter(b,1).XData);
    
    Ch4_Boxes.CenY(b) = mean(Ch4_gscatter(b,1).YData);
    Ch4_Boxes.CenX(b) = mean(Ch4_gscatter(b,1).XData);
    end
    
    Ch4_Cens(1:length(Ch4_Boxes.CenY),1) = Ch4_Boxes.CenY(:);
    Ch4_Cens(1:length(Ch4_Boxes.CenX),2) = Ch4_Boxes.CenX(:);
else end

%% Rangesearch Nearest Neighbor and compiling of "matched" VRC_IDs between channels.
clearvars VRC_RS RS_Idx RS_Idx2 RS_filt D;

if VRC_ID == [1,0,0,0],
    VRC_RS = table('Size',[Ch1_DBSCAN2_Max length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:length(Ch1_Boxes.MinY),
        VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(q);
        VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(q);
        VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(q);
        VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(q);
    end
elseif VRC_ID == [0,1,0,0],
    VRC_RS = table('Size',[Ch2_DBSCAN2_Max length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:length(Ch2_Boxes.MinY),
        VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(q);
        VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(q);
        VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(q);
        VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(q);
    end
elseif VRC_ID == [0,0,1,0],
    VRC_RS = table('Size',[Ch3_DBSCAN2_Max length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:length(Ch3_Boxes.MinY),
        VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(q);
        VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(q);
        VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(q);
        VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(q);
    end
elseif VRC_ID == [0,0,0,1],
    VRC_RS = table('Size',[Ch4_DBSCAN2_Max length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:length(Ch4_Boxes.MinY),
        VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(q);
        VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(q);
        VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(q);
        VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(q);
    end
elseif VRC_ID == [1,1,0,0],
   if VRC_IDCH == 1,
       [RS_Idx,D] = rangesearch(Ch2_Cens,Ch1_Cens,RS_Rad(1,1),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:size(RS_Idx,1),
         if RS_Idx{q,1} > 0,
        VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(RS_Idx{q,1});
        VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(RS_Idx{q,1});
        VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(RS_Idx{q,1});
        VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(RS_Idx{q,1});

        VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(q);
        VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(q);
        VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(q);
        VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(q);
         else end
         end
    VRC_RS = VRC_RS(RS_filt,:);
    
    elseif VRC_IDCH == 2,
        [RS_Idx,D] = rangesearch(Ch1_Cens,Ch2_Cens,RS_Rad(1,2),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    
            for q = 1:size(RS_Idx,1),
                if RS_Idx{q,1} > 0,
                VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(RS_Idx{q,1});
                VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(RS_Idx{q,1});
                VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(RS_Idx{q,1});
                VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(RS_Idx{q,1});

                VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(q);
                VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(q);
                VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(q);
                VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(q);
                else end
            end
            VRC_RS = VRC_RS(RS_filt,:);
        else end   

elseif VRC_ID == [1,0,1,0],
    if VRC_IDCH == 1,
       [RS_Idx,D] = rangesearch(Ch3_Cens,Ch1_Cens,RS_Rad(1,1),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:size(RS_Idx,1),
         if RS_Idx{q,1} > 0,
        VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(RS_Idx{q,1});
        VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(RS_Idx{q,1});
        VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(RS_Idx{q,1});
        VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(RS_Idx{q,1});

        VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(q);
        VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(q);
        VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(q);
        VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(q);
         else end
         end
    VRC_RS = VRC_RS(RS_filt,:);
    
    elseif VRC_IDCH == 3,
        [RS_Idx,D] = rangesearch(Ch1_Cens,Ch3_Cens,RS_Rad(1,3),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    
            for q = 1:size(RS_Idx,1),
                if RS_Idx{q,1} > 0,
                VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(RS_Idx{q,1});
                VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(RS_Idx{q,1});
                VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(RS_Idx{q,1});
                VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(RS_Idx{q,1});

                VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(q);
                VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(q);
                VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(q);
                VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(q);
                else end
            end
            VRC_RS = VRC_RS(RS_filt,:);
        else end   

elseif VRC_ID == [1,0,0,1],
   if VRC_IDCH == 1,
       [RS_Idx,D] = rangesearch(Ch4_Cens,Ch1_Cens,RS_Rad(1,1),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:size(RS_Idx,1),
         if RS_Idx{q,1} > 0,
        VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(RS_Idx{q,1});
        VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(RS_Idx{q,1});
        VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(RS_Idx{q,1});
        VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(RS_Idx{q,1});

        VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(q);
        VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(q);
        VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(q);
        VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(q);
         else end
         end
    VRC_RS = VRC_RS(RS_filt,:);
    
    elseif VRC_IDCH == 4,
        [RS_Idx,D] = rangesearch(Ch1_Cens,Ch4_Cens,RS_Rad(1,4),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    
            for q = 1:size(RS_Idx,1),
                if RS_Idx{q,1} > 0,
                VRC_RS.Ch1MinY(q) = Ch1_Boxes.MinY(RS_Idx{q,1});
                VRC_RS.Ch1MaxY(q) = Ch1_Boxes.MaxY(RS_Idx{q,1});
                VRC_RS.Ch1MinX(q) = Ch1_Boxes.MinX(RS_Idx{q,1});
                VRC_RS.Ch1MaxX(q) = Ch1_Boxes.MaxX(RS_Idx{q,1});

                VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(q);
                VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(q);
                VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(q);
                VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(q);
                else end
            end
            VRC_RS = VRC_RS(RS_filt,:);
        else end  
        
elseif VRC_ID == [0,1,1,0];
    if VRC_IDCH == 2,
       [RS_Idx,D] = rangesearch(Ch3_Cens,Ch2_Cens,RS_Rad(1,2),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:size(RS_Idx,1),
         if RS_Idx{q,1} > 0,
        VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(RS_Idx{q,1});
        VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(RS_Idx{q,1});
        VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(RS_Idx{q,1});
        VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(RS_Idx{q,1});

        VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(q);
        VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(q);
        VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(q);
        VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(q);
         else end
         end
    VRC_RS = VRC_RS(RS_filt,:);
    
    elseif VRC_IDCH == 3,
        [RS_Idx,D] = rangesearch(Ch2_Cens,Ch3_Cens,RS_Rad(1,3),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    
            for q = 1:size(RS_Idx,1),
                if RS_Idx{q,1} > 0,
                VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(RS_Idx{q,1});
                VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(RS_Idx{q,1});
                VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(RS_Idx{q,1});
                VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(RS_Idx{q,1});

                VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(q);
                VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(q);
                VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(q);
                VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(q);
                else end
            end
            VRC_RS = VRC_RS(RS_filt,:);
        else end   

elseif VRC_ID == [0,1,0,1];
    if VRC_IDCH == 2,
       [RS_Idx,D] = rangesearch(4,Ch2_Cens,RS_Rad(1,2),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:size(RS_Idx,1),
         if RS_Idx{q,1} > 0,
        VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(RS_Idx{q,1});
        VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(RS_Idx{q,1});
        VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(RS_Idx{q,1});
        VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(RS_Idx{q,1});

        VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(q);
        VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(q);
        VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(q);
        VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(q);
         else end
         end
    VRC_RS = VRC_RS(RS_filt,:);
    
    elseif VRC_IDCH == 4,
        [RS_Idx,D] = rangesearch(Ch2_Cens,Ch4_Cens,RS_Rad(1,4),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    
            for q = 1:size(RS_Idx,1),
                if RS_Idx{q,1} > 0,
                VRC_RS.Ch2MinY(q) = Ch2_Boxes.MinY(RS_Idx{q,1});
                VRC_RS.Ch2MaxY(q) = Ch2_Boxes.MaxY(RS_Idx{q,1});
                VRC_RS.Ch2MinX(q) = Ch2_Boxes.MinX(RS_Idx{q,1});
                VRC_RS.Ch2MaxX(q) = Ch2_Boxes.MaxX(RS_Idx{q,1});

                VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(q);
                VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(q);
                VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(q);
                VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(q);
                else end
            end
            VRC_RS = VRC_RS(RS_filt,:);
        else end   

elseif VRC_ID == [0,0,1,1];
   if VRC_IDCH == 3,
       [RS_Idx,D] = rangesearch(Ch4_Cens,Ch3_Cens,RS_Rad(1,3),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    for q = 1:size(RS_Idx,1),
         if RS_Idx{q,1} > 0,
        VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(RS_Idx{q,1});
        VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(RS_Idx{q,1});
        VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(RS_Idx{q,1});
        VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(RS_Idx{q,1});

        VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(q);
        VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(q);
        VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(q);
        VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(q);
         else end
         end
    VRC_RS = VRC_RS(RS_filt,:);
    
    elseif VRC_IDCH == 4,
        [RS_Idx,D] = rangesearch(Ch3_Cens,Ch4_Cens,RS_Rad(1,4),'SortIndices',true);
        
        for rs = 1:size(RS_Idx,1),
            if isempty(RS_Idx{rs,1}) == 1
                RS_Idx{rs,1} = 0;
            else end
            if length(RS_Idx{rs,1}) > 1,
                RS_Idx{rs,1} = RS_Idx{rs,1}(1,1);
            else end
        end
        RS_Idx2 = table2array(cell2table(RS_Idx));
        RS_filt = RS_Idx2>0;
        VRC_RS = table('Size',[size(RS_Idx,1) length(varNames_RS)],'VariableTypes',varTypes_RS,'VariableNames',varNames_RS);
    
            for q = 1:size(RS_Idx,1),
                if RS_Idx{q,1} > 0,
                VRC_RS.Ch3MinY(q) = Ch3_Boxes.MinY(RS_Idx{q,1});
                VRC_RS.Ch3MaxY(q) = Ch3_Boxes.MaxY(RS_Idx{q,1});
                VRC_RS.Ch3MinX(q) = Ch3_Boxes.MinX(RS_Idx{q,1});
                VRC_RS.Ch3MaxX(q) = Ch3_Boxes.MaxX(RS_Idx{q,1});

                VRC_RS.Ch4MinY(q) = Ch4_Boxes.MinY(q);
                VRC_RS.Ch4MaxY(q) = Ch4_Boxes.MaxY(q);
                VRC_RS.Ch4MinX(q) = Ch4_Boxes.MinX(q);
                VRC_RS.Ch4MaxX(q) = Ch4_Boxes.MaxX(q);
                else end
            end
            VRC_RS = VRC_RS(RS_filt,:);
        else end   

elseif sum(VRC_ID(:)) >2, disp(VRC_ID_high);
elseif sum(VRC_ID(:)) <1, disp(VRC_ID_low);
end

%% Creation of expanded crop boxes for each VRC, based on "matched" IDs in the previous section.
VRC_crops = table('Size',[size(VRC_RS,1) length(varNames_VRC)],'VariableTypes',varTypes_VRC,'VariableNames',varNames_VRC);
if sum(VRC_ID(:)) < 2,
 for v = 1:size(VRC_RS,1),
    VRC_crops.VRC_MinY(v) = table2array(VRC_RS(v,VRC_Idx2(1,1)));
    VRC_crops.VRC_MaxY(v) = table2array(VRC_RS(v,VRC_Idx2(2,1)));
    VRC_crops.VRC_MinX(v) = table2array(VRC_RS(v,VRC_Idx2(3,1)));
    VRC_crops.VRC_MaxX(v) = table2array(VRC_RS(v,VRC_Idx2(4,1)));
 end
else
    for v = 1:size(VRC_RS,1),
 VRC_crops.VRC_MinY(v) = min(table2array(VRC_RS(v,VRC_Idx2(1,1))),table2array(VRC_RS(v,VRC_Idx2(1,2))));
 VRC_crops.VRC_MaxY(v) = max(table2array(VRC_RS(v,VRC_Idx2(2,1))),table2array(VRC_RS(v,VRC_Idx2(2,2))));
 VRC_crops.VRC_MinX(v) = min(table2array(VRC_RS(v,VRC_Idx2(3,1))),table2array(VRC_RS(v,VRC_Idx2(3,2))));
 VRC_crops.VRC_MaxX(v) = max(table2array(VRC_RS(v,VRC_Idx2(4,1))),table2array(VRC_RS(v,VRC_Idx2(4,2))));
    end
end

%% Draw VRCs on image and save nuclear image as PDF
for slash = 1:length(strfind(filename,'\'))
    if slash ==1, filename2 = extractAfter(filename,'\');
    else filename2 = extractAfter(filename2,'\');
    end
end

filename2 = filename2(1:end-4); 
filename3 = append(filename2,'.pdf');
Folder_Outputs = append(Folder,'Outputs\');
    
if f==1, mkdir Outputs; cd(Folder_Outputs);
else cd(Folder_Outputs);
end
    
if drawVRCimages == 1,
    if VRC_ID == [1,0,0,0], imshow(Ch1(:,:,Ch1_MaxSlice(f,1)).*3);  %CHANGED from Ch1_seg to Ch1
    elseif VRC_ID == [0,1,0,0], imshow(Ch2_seg(:,:,Ch2_MaxSlice(f,1)).*3);
    elseif VRC_ID == [0,0,1,0], imshow(Ch3_seg(:,:,Ch3_MaxSlice(f,1)).*3);
    elseif VRC_ID == [0,0,0,1], imshow(Ch4_seg(:,:,Ch4_MaxSlice(f,1)).*3);
    elseif VRC_ID == [1,1,0,0], imshow(imfuse(Ch1_seg(:,:,Ch1_MaxSlice(f,1)).*3,Ch2_seg(:,:,Ch2_MaxSlice(f,1)).*3,'ColorChannels',[2 1 1]));
    elseif VRC_ID == [1,0,1,0], imshow(imfuse(Ch1_seg(:,:,Ch1_MaxSlice(f,1)).*3,Ch3_seg(:,:,Ch3_MaxSlice(f,1)).*3,'ColorChannels',[1 2 0]));
    elseif VRC_ID == [1,0,0,1], imshow(imfuse(Ch1_seg(:,:,Ch1_MaxSlice(f,1)).*3,Ch4_seg(:,:,Ch4_MaxSlice(f,1)).*3,'ColorChannels',[2 1 1]));
    elseif VRC_ID == [0,1,1,0], imshow(imfuse(Ch2_seg(:,:,Ch2_MaxSlice(f,1)).*3,Ch3_seg(:,:,Ch3_MaxSlice(f,1)).*3,'ColorChannels',[2 1 0]));
    elseif VRC_ID == [0,1,0,1], imshow(imfuse(Ch2_seg(:,:,Ch2_MaxSlice(f,1)).*3,Ch4_seg(:,:,Ch4_MaxSlice(f,1)).*3,'ColorChannels',[1 2 1]));
    elseif VRC_ID == [0,0,1,1], imshow(imfuse(Ch3_seg(:,:,Ch3_MaxSlice(f,1)).*3,Ch4_seg(:,:,Ch4_MaxSlice(f,1)).*3,'ColorChannels',[1 2 1]));
    else end
    
    hold on
    
    for a = 1:size(VRC_crops,1),
        rectangle('Position',[VRC_crops.VRC_MinX(a) VRC_crops.VRC_MinY(a) (VRC_crops.VRC_MaxX(a)-VRC_crops.VRC_MinX(a)) (VRC_crops.VRC_MaxY(a)-VRC_crops.VRC_MinY(a))],'EdgeColor','y','LineWidth',0.1); %[x y w h]
    end
    
    hold off
    ax = gca;
   
    exportgraphics(ax,filename3,'ContentType','vector','Resolution','500');

else
end

Settings = table({'DBSCAN_Rad';'DBSCAN_MinPts';'RS_Rad'},[DBSCAN_Rad(1,1);DBSCAN_MinPts(1,1);RS_Rad(1,1)],[DBSCAN_Rad(1,2);DBSCAN_MinPts(1,2);RS_Rad(1,2)],[DBSCAN_Rad(1,3);DBSCAN_MinPts(1,3);RS_Rad(1,3)],[DBSCAN_Rad(1,4);DBSCAN_MinPts(1,4);RS_Rad(1,4)],'VariableNames',{'Variable','Ch1','Ch2','Ch3','Ch4'});
cd(Folder_Outputs), writetable(Settings,'DBSCAN-VariableValues.xlsx','WriteVariableNames',true); cd(Folder);
%% Dissect individual VRCs and save images (if drawVRCimages == 1).

    cd(Folder_Outputs);
    Folder_VRCs=append(filename2,' VRCs\');

if drawVRCimages == 1,
    for t = 1:size(VRC_crops,1),
        if t == 1, mkdir(Folder_VRCs); cd(Folder_VRCs);
            delete *.tif;
        else end
        
        VRC_subname = num2str(t);
        if t <10, VRC_filename = append(filename2,' VRC 0',VRC_subname,'.tiff');
        else VRC_filename = append(filename2,' VRC ',VRC_subname,'.tiff');
        end
        
        VRCs{t,1}(:,:,:,1) = Ch1(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
        VRCs{t,1}(:,:,:,2) = Ch2(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
        if Channels >2, VRCs{t,1}(:,:,:,3) = Ch3_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
        if Channels >3, VRCs{t,1}(:,:,:,4) = Ch4_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
                
        bfsave(VRCs{t,1}(:,:,:,:),VRC_filename,'dimensionOrder','XYZTC');
    end
    
else 
    for t = 1:size(VRC_crops,1),       
        VRCs{t,1}(:,:,:,1) = Ch1_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
        VRCs{t,1}(:,:,:,2) = Ch2_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:);
        if Channels >2, VRCs{t,1}(:,:,:,3) = Ch3_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
        if Channels >3, VRCs{t,1}(:,:,:,4) = Ch4_seg(VRC_crops.VRC_MinY(t):VRC_crops.VRC_MaxY(t), VRC_crops.VRC_MinX(t):VRC_crops.VRC_MaxX(t),:); else end;
    end
 
end

%% Volume Analysis for each VRC %%
% VRCs{#,1}(X,Y,Z,C).

if AnalysisVolume == 1
    
clearvars VRC_Volumes;
VRC_Volumes = table('Size',[size(VRCs,1) length(varNames_Volumes)],'VariableTypes',varTypes_Volumes,'VariableNames',varNames_Volumes);

cd(Folder_Outputs);
if f == 1, mkdir('Volume Analysis'); cd('Volume Analysis');
    else cd('Volume Analysis'); 
end

for h = 1:size(VRCs,1),
    
    Volumes_filename = append(filename2,' Volumes.xlsx');

    for u = 1:size(VRCs,1),
        VRC_Volumes.Ch1Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,1)>0)));
        VRC_Volumes.Ch2Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,2)>0)));
        if Channels>2, VRC_Volumes.Ch3Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,3)>0))); else, VRC_Volumes.Ch3Pix(:) = 0; end
        if Channels>3, VRC_Volumes.Ch4Pix(u) = sum(sum(sum(VRCs{u,1}(:,:,:,4)>0))); else, VRC_Volumes.Ch4Pix(:) = 0; end
    end
end

writetable(VRC_Volumes,Volumes_filename,'WriteVariableNames',true); cd(Folder);

else end
%% PCC Analysis for each VRC %%

if AnalysisPCC ==1,

clearvars VRC_PCCs;
VRC_PCCs = table('Size',[size(VRCs,1) length(varNames_PCCs)],'VariableTypes',varTypes_PCCs,'VariableNames',varNames_PCCs);

cd(Folder_Outputs);
if f == 1, mkdir('PCC Analysis'); cd('PCC Analysis');
    else cd('PCC Analysis');
end

for hh = 1:size(VRCs,1),
    
    PCCs_filename = append(filename2,' PCCs.xlsx');
    VRC_area = (VRC_crops.VRC_MaxY(1) - VRC_crops.VRC_MinY(1))*((VRC_crops.MaxX(1) - VRC_crops.VRC_minX(1)));
    
    for uu = 1:size(VRCs,1),
        Ch1_re = reshape(Ch1_seg,[VRC_area 1]);
        
    end
end
else end
end
close all %Closes all open figures that are made from "gscatter". Not sure if necessarily the best way, but it works.
toc




