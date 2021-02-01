function [Percentile_Values,BWMask,Seg] = PercentileThresholding(f,I,Channels,Percentile_Thresholds)

%% Determine the intensity values of the percentile ranges provided in Percentile_Thresholds %%
Percentile_Values.Lower(1,1) = prctile(I.Ch1,Percentile_Thresholds.Lower(1,1),'all');
Percentile_Values.Upper(1,1) = prctile(I.Ch1,Percentile_Thresholds.Upper(1,1),'all');
Percentile_Values.Lower(1,2) = prctile(I.Ch2,Percentile_Thresholds.Lower(1,2),'all');
Percentile_Values.Upper(1,2) = prctile(I.Ch2,Percentile_Thresholds.Upper(1,2),'all');
if Channels > 2
    Percentile_Values.Lower(1,3) = prctile(I.Ch3,Percentile_Thresholds.Lower(1,3),'all');
    Percentile_Values.Upper(1,3) = prctile(I.Ch3,Percentile_Thresholds.Upper(1,3),'all');
else end
if Channels > 3
    Percentile_Values.Lower(1,4) = prctile(I.Ch4,Percentile_Thresholds.Lower(1,4),'all');
    Percentile_Values.Upper(1,4) = prctile(I.Ch4,Percentile_Thresholds.Upper(1,4),'all');
else end

%% Mask and Segment each channel based on percentile values %%

BWMask.Ch1 = (Percentile_Values.Upper(1,1)>=I.Ch1) & (I.Ch1>Percentile_Values.Lower(1,1));
Seg.Ch1 = uint16(double(I.Ch1).*BWMask.Ch1);

BWMask.Ch2 = (Percentile_Values.Upper(1,2)>=I.Ch2) & (I.Ch2>Percentile_Values.Lower(1,2));
Seg.Ch2 = uint16(double(I.Ch2).*BWMask.Ch2);

if Channels > 2
    BWMask.Ch3 = (Percentile_Values.Upper(1,3)>=I.Ch3) & (I.Ch3>Percentile_Values.Lower(1,3));
    Seg.Ch3 = uint16(double(I.Ch3).*BWMask.Ch3);
else end

if Channels > 3
    BWMask.Ch4 = (Percentile_Values.Upper(1,4)>=I.Ch4) & (I.Ch4>Percentile_Values.Lower(1,4));
    Seg.Ch4 = uint16(double(I.Ch4).*BWMask.Ch4);
else end

end

