
function [Results_NuclearAnalysis] = NuclearAnalysis(DAPI_Watershed_BW2,Channels,DAPI,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4)
 
Results_NuclearAnalysis.NuclearImage = DAPI_Watershed_BW2;
Results_NuclearAnalysis.NucleiCC = bwconncomp(Results_NuclearAnalysis.NuclearImage,8);
Results_NuclearAnalysis.NuclearProps = regionprops(Results_NuclearAnalysis.NucleiCC,DAPI,'Area','Circularity','Centroid');

Results_NuclearAnalysis.NucleiNumber = length(Results_NuclearAnalysis.NuclearProps);

for p = 1:Results_NuclearAnalysis.NucleiCC.NumObjects
    Results_NuclearAnalysis.MeanCh1(p,1) = mean(Focus_Ch1(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
    Results_NuclearAnalysis.SumCh1(p,1) = sum(Focus_Ch1(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
    Results_NuclearAnalysis.MeanCh2(p,1) = mean(Focus_Ch2(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
    Results_NuclearAnalysis.SumCh2(p,1) = sum(Focus_Ch2(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
    if Channels>2
        Results_NuclearAnalysis.MeanCh3(p,1) = mean(Focus_Ch3(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
        Results_NuclearAnalysis.SumCh3(p,1) = sum(Focus_Ch3(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
    else
    end
    if Channels>3
        Results_NuclearAnalysis.MeanCh4(p,1) = mean(Focus_Ch4(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
        Results_NuclearAnalysis.SumCh4(p,1) = sum(Focus_Ch4(Results_NuclearAnalysis.NucleiCC.PixelIdxList{p}));
    else
    end
end

