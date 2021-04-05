
function [DAPI_Watershed_BW2,DAPI_Watershed_Perim,DAPIseg_cc,DAPIseg_props] = NucleiSegmentation(Ch_DAPI,ResY)

I_DAPI_bg = imopen(Ch_DAPI,strel('disk',(round(ResY/50))));
I_DAPI2 = Ch_DAPI - I_DAPI_bg;
I_DAPI2_BW = imfill(imbinarize(imadjust(I_DAPI2)),'holes');

DAPI_DistTransf1 = bwdist(~I_DAPI2_BW);
DAPI_DistTransf1 = imcomplement(DAPI_DistTransf1);
DAPI_DistTransfMask = imextendedmin(DAPI_DistTransf1,2);
DAPI_DistTransf2 = imimposemin(DAPI_DistTransf1,DAPI_DistTransfMask);
DAPI_Watershed = watershed(DAPI_DistTransf2);
DAPI_Watershed_BW = I_DAPI2_BW;
DAPI_Watershed_BW(DAPI_Watershed == 0) = 0;

DAPI_Watershed_BW2 = logical(imdilate(imerode(DAPI_Watershed_BW,strel('disk',1)),strel('disk',1)));
DAPI_Watershed_BW2 = bwareafilt(DAPI_Watershed_BW2,[1000 10000]);

DAPI_Watershed_Perim = bwperim(imdilate(DAPI_Watershed_BW2,strel('disk',1)));

DAPIseg_cc = bwconncomp(DAPI_Watershed_BW2,8);
DAPIseg_props = regionprops(DAPIseg_cc,Ch_DAPI,'MeanIntensity','PixelValues','Area','Centroid','Circularity');

for p = 1:length(DAPIseg_props)
DAPI_NuclearAreas(p,1) = DAPIseg_props(p).Area;
end
end

