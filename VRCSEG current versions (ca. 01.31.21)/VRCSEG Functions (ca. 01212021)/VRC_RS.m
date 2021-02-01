function [RSInfo] = VRC_RS(f,RS_Rad,VRC_ID,VRC_RSseed,ClusterInfo,RSInfo)
%% RangeSearch Analysis %%

    if VRC_ID(1,1) == 1 && VRC_RSseed == 1
        RSInfo(f).SeedBoundingBoxes = ClusterInfo(f).Ch1.BoundingBoxes;
        RSInfo(f).SeedCentroidCoords = [ClusterInfo(f).Ch1.ClusterCentroidsY ClusterInfo(f).Ch1.ClusterCentroidsX];
        RSInfo(f).SeedClusterNumber = size(ClusterInfo(f).Ch1.BoundingBoxes,1);
        RSInfo(f).SeedChannel = 1;
    elseif VRC_ID(1,1) == 1 && VRC_RSseed ~=1
        RSInfo(f).SearchBoundingBoxes = ClusterInfo(f).Ch1.BoundingBoxes;
        RSInfo(f).SearchCentroidCoords = [ClusterInfo(f).Ch1.ClusterCentroidsY ClusterInfo(f).Ch1.ClusterCentroidsX];
        RSInfo(f).SearchChannel = 1;
    else end
    
    if VRC_ID(1,2) == 1 && VRC_RSseed == 2
        RSInfo(f).SeedBoundingBoxes = ClusterInfo(f).Ch2.BoundingBoxes;
        RSInfo(f).SeedCentroidCoords = [ClusterInfo(f).Ch2.ClusterCentroidsY ClusterInfo(f).Ch2.ClusterCentroidsX];
        RSInfo(f).SeedClusterNumber = size(ClusterInfo(f).Ch2.BoundingBoxes,1);
        RSInfo(f).SeedChannel = 2;
    elseif VRC_ID(1,2) == 1 && VRC_RSseed ~=2
        RSInfo(f).SearchBoundingBoxes = ClusterInfo(f).Ch2.BoundingBoxes;
        RSInfo(f).SearchCentroidCoords = [ClusterInfo(f).Ch2.ClusterCentroidsY ClusterInfo(f).Ch2.ClusterCentroidsX];
        RSInfo(f).SearchChannel = 2;
    else end
     
    if VRC_ID(1,3) == 1 && VRC_RSseed == 3
        RSInfo(f).SeedBoundingBoxes = ClusterInfo(f).Ch3.BoundingBoxes;
        RSInfo(f).SeedCentroidCoords = [ClusterInfo(f).Ch3.ClusterCentroidsY ClusterInfo(f).Ch3.ClusterCentroidsX];
        RSInfo(f).SeedClusterNumber = size(ClusterInfo(f).Ch3.BoundingBoxes,1);
        RSInfo(f).SeedChannel = 3;
    elseif VRC_ID(1,3) == 1 && VRC_RSseed ~=3
        RSInfo(f).SearchBoundingBoxes = ClusterInfo(f).Ch3.BoundingBoxes;
        RSInfo(f).SearchCentroidCoords = [ClusterInfo(f).Ch3.ClusterCentroidsY ClusterInfo(f).Ch3.ClusterCentroidsX];
        RSInfo(f).SearchChannel = 3;
    else end
    
    if VRC_ID(1,4) == 1 && VRC_RSseed == 4
        RSInfo(f).SeedBoundingBoxes = ClusterInfo(f).Ch4.BoundingBoxes;
        RSInfo(f).SeedCentroidCoords = [ClusterInfo(f).Ch4.ClusterCentroidsY ClusterInfo(f).Ch4.ClusterCentroidsX];
        RSInfo(f).SeedClusterNumber = size(ClusterInfo(f).Ch4.BoundingBoxes,1);
        RSInfo(f).SeedChannel = 4;
    elseif VRC_ID(1,4) == 1 && VRC_RSseed ~=4
        RSInfo(f).SearchBoundingBoxes = ClusterInfo(f).Ch4.BoundingBoxes;
        RSInfo(f).SearchCentroidCoords = [ClusterInfo(f).Ch4.ClusterCentroidsY ClusterInfo(f).Ch4.ClusterCentroidsX];
        RSInfo(f).SearchChannel = 4;
    else end
    
    RSInfo(f).SearchRadius = RS_Rad(RSInfo(f).SeedChannel);
    
    if sum(VRC_ID,'all') == 2
    [RSInfo(f).Idx,RSInfo(f).Distance] = rangesearch(RSInfo(f).SearchCentroidCoords,RSInfo(f).SeedCentroidCoords,RSInfo(f).SearchRadius,'SortIndices',true);
    for s = 1:size(RSInfo(f).Idx,1)
        if isempty(RSInfo(f).Idx{s,1}) == 1
            RSInfo(f).Idx{s,1} = 0;
        else end
        if length(RSInfo(f).Idx{s,1}) > 1
            RSInfo(f).Idx{s,1} = RSInfo(f).Idx{s,1}(1,1);
        else end
    end
    RSInfo(f).IdxArray = table2array(cell2table(RSInfo(f).Idx));
    RSInfo(f).IdxFilter = RSInfo(f).IdxArray>0;
    else end
    
    for q = 1:size(RSInfo(f).SeedCentroidCoords,1)
        if sum(VRC_ID,'all') == 2
            if RSInfo(f).Idx{q,1} > 0
            RSInfo(f).SearchCoords(q).MinY = RSInfo(f).SearchBoundingBoxes(RSInfo(f).Idx{q,1}).MinY;
            RSInfo(f).SearchCoords(q).MaxY = RSInfo(f).SearchBoundingBoxes(RSInfo(f).Idx{q,1}).MaxY;
            RSInfo(f).SearchCoords(q).MinX = RSInfo(f).SearchBoundingBoxes(RSInfo(f).Idx{q,1}).MinX;
            RSInfo(f).SearchCoords(q).MaxX = RSInfo(f).SearchBoundingBoxes(RSInfo(f).Idx{q,1}).MaxX;
            else end
        else end
        RSInfo(f).SeedCoords(q).MinY = RSInfo(f).SeedBoundingBoxes(q).MinY;
        RSInfo(f).SeedCoords(q).MaxY = RSInfo(f).SeedBoundingBoxes(q).MaxY;
        RSInfo(f).SeedCoords(q).MinX = RSInfo(f).SeedBoundingBoxes(q).MinX;
        RSInfo(f).SeedCoords(q).MaxX = RSInfo(f).SeedBoundingBoxes(q).MaxX;
    end
    
    if sum(VRC_ID,'all') == 2
    RSInfo(f).SeedCoords = RSInfo(f).SeedCoords(RSInfo(f).IdxFilter)';
    RSInfo(f).SearchCoords = RSInfo(f).SearchCoords(RSInfo(f).IdxFilter)';
    else
    RSInfo(f).SeedCoords = RSInfo(f).SeedCoords';
    end
 
%% Create "expanded" crop boxes for each VRC, based on RangeSearch results. %%
if sum(VRC_ID,'all') == 1
    for v = 1:size(RSInfo(f).SeedCoords)
        RSInfo(f).ExpandedBB(v,1).MinY = RSInfo(f).SeedCoords(v,1).MinY;
        RSInfo(f).ExpandedBB(v,1).MaxY = RSInfo(f).SeedCoords(v,1).MaxY;
        RSInfo(f).ExpandedBB(v,1).MinX = RSInfo(f).SeedCoords(v,1).MinX;
        RSInfo(f).ExpandedBB(v,1).MaxX = RSInfo(f).SeedCoords(v,1).MaxX;
    end
elseif sum(VRC_ID,'all') == 2
    for v = 1:size(RSInfo(f).SeedCoords)
        RSInfo(f).ExpandedBB(v,1).MinY = min(RSInfo(f).SeedCoords(v,1).MinY,RSInfo(f).SearchCoords(v,1).MinY);
        RSInfo(f).ExpandedBB(v,1).MaxY = max(RSInfo(f).SeedCoords(v,1).MaxY,RSInfo(f).SearchCoords(v,1).MaxY);
        RSInfo(f).ExpandedBB(v,1).MinX = min(RSInfo(f).SeedCoords(v,1).MinX,RSInfo(f).SearchCoords(v,1).MinX);
        RSInfo(f).ExpandedBB(v,1).MaxX = max(RSInfo(f).SeedCoords(v,1).MaxX,RSInfo(f).SearchCoords(v,1).MaxX);
    end
else
end
    
end