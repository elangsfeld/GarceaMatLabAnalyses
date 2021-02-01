function [ClusterInfo,DBSCAN] = VRC_DBSCAN(f,DBSCAN_Rad,DBSCAN_MinPts,Channels,VRC_ID,FocusInfo,BWMask,ClusterInfo,DBSCAN)

%% If Channel #1 is used for identifying VRCs. %%

if VRC_ID(1,1) == 1
    [ClusterInfo(f).Ch1.Row,ClusterInfo(f).Ch1.Column] = find(BWMask.Ch1(:,:,FocusInfo(f).MaxPlane(1,1))); %IDs above-threshold pixsels in the brightest z-slice.
    ClusterInfo(f).Ch1.XYCoord = [ClusterInfo(f).Ch1.Column, ClusterInfo(f).Ch1.Row]; %Combines positive coordinates into one array.
    DBSCAN(f).Ch1.Idx = dbscan(ClusterInfo(f).Ch1.XYCoord,DBSCAN_Rad(1,1),DBSCAN_MinPts(1,1)); %Performs DBSCAN analysis based on input values at top of script.
    DBSCAN(f).Ch1.IdxFilter = DBSCAN(f).Ch1.Idx>0; %Makes logical for DBSCAN groups (>0 ID value). Outliers get -1 ID value that needs to be filtered out.
    DBSCAN(f).Ch1.Filtered = DBSCAN(f).Ch1.Idx(DBSCAN(f).Ch1.IdxFilter); %Filters original DBSCAN to remove -1 ID values.
    DBSCAN(f).Ch1.ClusterNumber = max(DBSCAN(f).Ch1.Filtered(:)); %Determines how many clusters were detected.
    ClusterInfo(f).Ch1.RowFiltered = ClusterInfo(f).Ch1.Row(DBSCAN(f).Ch1.IdxFilter); %Filters rows to remove -1 ID values.
    ClusterInfo(f).Ch1.ColumnFiltered = ClusterInfo(f).Ch1.Column(DBSCAN(f).Ch1.IdxFilter); %Filters columns to remove -1 ID values.
    ClusterInfo(f).Ch1.XYCoordFiltered = [ClusterInfo(f).Ch1.ColumnFiltered, ClusterInfo(f).Ch1.RowFiltered]; %Places filtered Row/Col into array.
    
    ClusterInfo(f).Ch1.gscatter = gscatter(ClusterInfo(f).Ch1.XYCoordFiltered(:,1),ClusterInfo(f).Ch1.XYCoordFiltered(:,2),DBSCAN(f).Ch1.Filtered,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.
    
    for b = 1:DBSCAN(f).Ch1.ClusterNumber %This loop finds the min/max Y/X values for each "cluster IDed by gscatter fxn above, then calculates their centroid(s).
        ClusterInfo(f).Ch1.BoundingBoxes(b,1).MinY = min(ClusterInfo(f).Ch1.gscatter(b,1).YData);
        ClusterInfo(f).Ch1.BoundingBoxes(b,1).MaxY = max(ClusterInfo(f).Ch1.gscatter(b,1).YData);
        ClusterInfo(f).Ch1.BoundingBoxes(b,1).MinX = min(ClusterInfo(f).Ch1.gscatter(b,1).XData);
        ClusterInfo(f).Ch1.BoundingBoxes(b,1).MaxX = max(ClusterInfo(f).Ch1.gscatter(b,1).XData);
        ClusterInfo(f).Ch1.BoundingBoxes(b,1).CentroidY = mean(ClusterInfo(f).Ch1.gscatter(b,1).YData);
        ClusterInfo(f).Ch1.ClusterCentroidsY(b,1) = ClusterInfo(f).Ch1.BoundingBoxes(b,1).CentroidY;
        ClusterInfo(f).Ch1.BoundingBoxes(b,1).CentroidX = mean(ClusterInfo(f).Ch1.gscatter(b,1).XData);
        ClusterInfo(f).Ch1.ClusterCentroidsX(b,1) = ClusterInfo(f).Ch1.BoundingBoxes(b,1).CentroidX;
    end    
else
end

%% If Channel #2 is used for identifying VRCs. %%

if VRC_ID(1,2) == 1
    [ClusterInfo(f).Ch2.Row,ClusterInfo(f).Ch2.Column] = find(BWMask.Ch2(:,:,FocusInfo(f).MaxPlane(1,1))); %IDs above-threshold pixsels in the brightest z-slice.
    ClusterInfo(f).Ch2.XYCoord = [ClusterInfo(f).Ch2.Column, ClusterInfo(f).Ch2.Row]; %Combines positive coordinates into one array.
    DBSCAN(f).Ch2.Idx = dbscan(ClusterInfo(f).Ch2.XYCoord,DBSCAN_Rad(1,1),DBSCAN_MinPts(1,1)); %Performs DBSCAN analysis based on input values at top of script.
    DBSCAN(f).Ch2.IdxFilter = DBSCAN(f).Ch2.Idx>0; %Makes logical for DBSCAN groups (>0 ID value). Outliers get -1 ID value that needs to be filtered out.
    DBSCAN(f).Ch2.Filtered = DBSCAN(f).Ch2.Idx(DBSCAN(f).Ch2.IdxFilter); %Filters original DBSCAN to remove -1 ID values.
    DBSCAN(f).Ch2.ClusterNumber = max(DBSCAN(f).Ch2.Filtered(:)); %Determines how many clusters were detected.
    ClusterInfo(f).Ch2.RowFiltered = ClusterInfo(f).Ch2.Row(DBSCAN(f).Ch2.IdxFilter); %Filters rows to remove -1 ID values.
    ClusterInfo(f).Ch2.ColumnFiltered = ClusterInfo(f).Ch2.Column(DBSCAN(f).Ch2.IdxFilter); %Filters columns to remove -1 ID values.
    ClusterInfo(f).Ch2.XYCoordFiltered = [ClusterInfo(f).Ch2.ColumnFiltered, ClusterInfo(f).Ch2.RowFiltered]; %Places filtered Row/Col into array.
    
    ClusterInfo(f).Ch2.gscatter = gscatter(ClusterInfo(f).Ch2.XYCoordFiltered(:,1),ClusterInfo(f).Ch2.XYCoordFiltered(:,2),DBSCAN(f).Ch2.Filtered,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.
    
    for b = 1:DBSCAN(f).Ch2.ClusterNumber %This loop finds the min/max Y/X values for each "cluster IDed by gscatter fxn above, then calculates their centroid(s).
        ClusterInfo(f).Ch2.BoundingBoxes(b,1).MinY = min(ClusterInfo(f).Ch2.gscatter(b,1).YData);
        ClusterInfo(f).Ch2.BoundingBoxes(b,1).MaxY = max(ClusterInfo(f).Ch2.gscatter(b,1).YData);
        ClusterInfo(f).Ch2.BoundingBoxes(b,1).MinX = min(ClusterInfo(f).Ch2.gscatter(b,1).XData);
        ClusterInfo(f).Ch2.BoundingBoxes(b,1).MaxX = max(ClusterInfo(f).Ch2.gscatter(b,1).XData);
        ClusterInfo(f).Ch2.BoundingBoxes(b,1).CentroidY = mean(ClusterInfo(f).Ch2.gscatter(b,1).YData);
        ClusterInfo(f).Ch2.ClusterCentroidsY(b,1) = ClusterInfo(f).Ch2.BoundingBoxes(b,1).CentroidY;
        ClusterInfo(f).Ch2.BoundingBoxes(b,1).CentroidX = mean(ClusterInfo(f).Ch2.gscatter(b,1).XData);
        ClusterInfo(f).Ch2.ClusterCentroidsX(b,1) = ClusterInfo(f).Ch2.BoundingBoxes(b,1).CentroidX;
    end    
else
end

%% If Channel #3 is used for identifying VRCs. %%

if Channels>2 && VRC_ID(1,3) == 1
    [ClusterInfo(f).Ch3.Row,ClusterInfo(f).Ch3.Column] = find(BWMask.Ch3(:,:,FocusInfo(f).MaxPlane(1,1))); %IDs above-threshold pixsels in the brightest z-slice.
    ClusterInfo(f).Ch3.XYCoord = [ClusterInfo(f).Ch3.Column, ClusterInfo(f).Ch3.Row]; %Combines positive coordinates into one array.
    DBSCAN(f).Ch3.Idx = dbscan(ClusterInfo(f).Ch3.XYCoord,DBSCAN_Rad(1,1),DBSCAN_MinPts(1,1)); %Performs DBSCAN analysis based on input values at top of script.
    DBSCAN(f).Ch3.IdxFilter = DBSCAN(f).Ch3.Idx>0; %Makes logical for DBSCAN groups (>0 ID value). Outliers get -1 ID value that needs to be filtered out.
    DBSCAN(f).Ch3.Filtered = DBSCAN(f).Ch3.Idx(DBSCAN(f).Ch3.IdxFilter); %Filters original DBSCAN to remove -1 ID values.
    DBSCAN(f).Ch3.ClusterNumber = max(DBSCAN(f).Ch3.Filtered(:)); %Determines how many clusters were detected.
    ClusterInfo(f).Ch3.RowFiltered = ClusterInfo(f).Ch3.Row(DBSCAN(f).Ch3.IdxFilter); %Filters rows to remove -1 ID values.
    ClusterInfo(f).Ch3.ColumnFiltered = ClusterInfo(f).Ch3.Column(DBSCAN(f).Ch3.IdxFilter); %Filters columns to remove -1 ID values.
    ClusterInfo(f).Ch3.XYCoordFiltered = [ClusterInfo(f).Ch3.ColumnFiltered, ClusterInfo(f).Ch3.RowFiltered]; %Places filtered Row/Col into array.
    
    ClusterInfo(f).Ch3.gscatter = gscatter(ClusterInfo(f).Ch3.XYCoordFiltered(:,1),ClusterInfo(f).Ch3.XYCoordFiltered(:,2),DBSCAN(f).Ch3.Filtered,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.
    
    for b = 1:DBSCAN(f).Ch3.ClusterNumber %This loop finds the min/max Y/X values for each "cluster IDed by gscatter fxn above, then calculates their centroid(s).
        ClusterInfo(f).Ch3.BoundingBoxes(b,1).MinY = min(ClusterInfo(f).Ch3.gscatter(b,1).YData);
        ClusterInfo(f).Ch3.BoundingBoxes(b,1).MaxY = max(ClusterInfo(f).Ch3.gscatter(b,1).YData);
        ClusterInfo(f).Ch3.BoundingBoxes(b,1).MinX = min(ClusterInfo(f).Ch3.gscatter(b,1).XData);
        ClusterInfo(f).Ch3.BoundingBoxes(b,1).MaxX = max(ClusterInfo(f).Ch3.gscatter(b,1).XData);
        ClusterInfo(f).Ch3.BoundingBoxes(b,1).CentroidY = mean(ClusterInfo(f).Ch3.gscatter(b,1).YData);
        ClusterInfo(f).Ch3.ClusterCentroidsY(b,1) = ClusterInfo(f).Ch3.BoundingBoxes(b,1).CentroidY;
        ClusterInfo(f).Ch3.BoundingBoxes(b,1).CentroidX = mean(ClusterInfo(f).Ch3.gscatter(b,1).XData);
        ClusterInfo(f).Ch3.ClusterCentroidsX(b,1) = ClusterInfo(f).Ch3.BoundingBoxes(b,1).CentroidX;
    end    
elseif Channels<3
    warning('There must be at least 3 channels for Channel 3 to be used for VRC Identification by DBSCAN.');
    warning('Check the values of "Channels" and "VRC_ID" variables.');
else
end

%% If Channel #4 is used for identifying VRCs. %%

if Channels>3 && VRC_ID(1,4) == 1
    [ClusterInfo(f).Ch4.Row,ClusterInfo(f).Ch4.Column] = find(BWMask.Ch4(:,:,FocusInfo(f).MaxPlane(1,1))); %IDs above-threshold pixsels in the brightest z-slice.
    ClusterInfo(f).Ch4.XYCoord = [ClusterInfo(f).Ch4.Column, ClusterInfo(f).Ch4.Row]; %Combines positive coordinates into one array.
    DBSCAN(f).Ch4.Idx = dbscan(ClusterInfo(f).Ch4.XYCoord,DBSCAN_Rad(1,1),DBSCAN_MinPts(1,1)); %Performs DBSCAN analysis based on input values at top of script.
    DBSCAN(f).Ch4.IdxFilter = DBSCAN(f).Ch4.Idx>0; %Makes logical for DBSCAN groups (>0 ID value). Outliers get -1 ID value that needs to be filtered out.
    DBSCAN(f).Ch4.Filtered = DBSCAN(f).Ch4.Idx(DBSCAN(f).Ch4.IdxFilter); %Filters original DBSCAN to remove -1 ID values.
    DBSCAN(f).Ch4.ClusterNumber = max(DBSCAN(f).Ch4.Filtered(:)); %Determines how many clusters were detected.
    ClusterInfo(f).Ch4.RowFiltered = ClusterInfo(f).Ch4.Row(DBSCAN(f).Ch4.IdxFilter); %Filters rows to remove -1 ID values.
    ClusterInfo(f).Ch4.ColumnFiltered = ClusterInfo(f).Ch4.Column(DBSCAN(f).Ch4.IdxFilter); %Filters columns to remove -1 ID values.
    ClusterInfo(f).Ch4.XYCoordFiltered = [ClusterInfo(f).Ch4.ColumnFiltered, ClusterInfo(f).Ch4.RowFiltered]; %Places filtered Row/Col into array.
    
    ClusterInfo(f).Ch4.gscatter = gscatter(ClusterInfo(f).Ch4.XYCoordFiltered(:,1),ClusterInfo(f).Ch4.XYCoordFiltered(:,2),DBSCAN(f).Ch4.Filtered,[],[],[1],'off tight'); %Creates Line variable with each DBSCAN-IDed group.
    
    for b = 1:DBSCAN(f).Ch4.ClusterNumber %This loop finds the min/max Y/X values for each "cluster IDed by gscatter fxn above, then calculates their centroid(s).
        ClusterInfo(f).Ch4.BoundingBoxes(b,1).MinY = min(ClusterInfo(f).Ch4.gscatter(b,1).YData);
        ClusterInfo(f).Ch4.BoundingBoxes(b,1).MaxY = max(ClusterInfo(f).Ch4.gscatter(b,1).YData);
        ClusterInfo(f).Ch4.BoundingBoxes(b,1).MinX = min(ClusterInfo(f).Ch4.gscatter(b,1).XData);
        ClusterInfo(f).Ch4.BoundingBoxes(b,1).MaxX = max(ClusterInfo(f).Ch4.gscatter(b,1).XData);
        ClusterInfo(f).Ch4.BoundingBoxes(b,1).CentroidY = mean(ClusterInfo(f).Ch4.gscatter(b,1).YData);
        ClusterInfo(f).Ch4.ClusterCentroidsY(b,1) = ClusterInfo(f).Ch4.BoundingBoxes(b,1).CentroidY;
        ClusterInfo(f).Ch4.BoundingBoxes(b,1).CentroidX = mean(ClusterInfo(f).Ch4.gscatter(b,1).XData);
        ClusterInfo(f).Ch4.ClusterCentroidsX(b,1) = ClusterInfo(f).Ch4.BoundingBoxes(b,1).CentroidX;
    end    
elseif Channels<4 && VRC_ID(1,4) == 1
    warning('There must be 4 channels for Channel 4 to be used for VRC Identification by DBSCAN.');
    warning('Check the values of "Channels" and "VRC_ID" variables.');
else
end   
    
end