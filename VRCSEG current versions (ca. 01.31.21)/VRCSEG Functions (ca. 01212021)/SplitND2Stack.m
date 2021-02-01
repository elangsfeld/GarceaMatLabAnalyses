function [I,FocusInfo] = SplitND2Stack(Channels,Res,RawImage,f,FocusInfo)

Slices = (length(RawImage{1,1})/Channels);

for i = 1:Slices
    ChannelPlanes(i).Ch1 = 1+(Channels*i-Channels); 
    I.Ch1(:,:,i) = RawImage{1,1}{ChannelPlanes(i).Ch1,1};
    FocusInfo(f).IntDenPerPlane(i,1) = sum(sum(I.Ch1(:,:,i)));
    
    ChannelPlanes(i).Ch2 = 2+(Channels*i-Channels); 
    I.Ch2(:,:,i) = RawImage{1,1}{ChannelPlanes(i).Ch2,1};
    FocusInfo(f).IntDenPerPlane(i,2) = sum(sum(I.Ch2(:,:,i)));
    
    if Channels>2 
        ChannelPlanes(i).Ch3 = 3+(Channels*i-Channels);
        I.Ch3(:,:,i) = RawImage{1,1}{ChannelPlanes(i).Ch3,1};
        FocusInfo(f).IntDenPerPlane(i,3) = sum(sum(I.Ch3(:,:,i)));
    else
        I.Ch3 = zeros(Res,Res,Slices);
    end
    if Channels>3
        ChannelPlanes(i).Ch4 = 4+(Channels*i-Channels);
        I.Ch4(:,:,i) = RawImage{1,1}{ChannelPlanes(i).Ch4,1};
        FocusInfo(f).IntDenPerPlane(i,4) = sum(sum(I.Ch4(:,:,i)));
    else
        I.Ch4 = zeros(Res,Res,Slices);
    end
end

[FocusInfo(f).MaxValue(1,1) FocusInfo(f).MaxPlane(1,1)] = max(FocusInfo(f).IntDenPerPlane(:,1));
[FocusInfo(f).MaxValue(1,2) FocusInfo(f).MaxPlane(1,2)] = max(FocusInfo(f).IntDenPerPlane(:,2));

if Channels>2
    [FocusInfo(f).MaxValue(1,3) FocusInfo(f).MaxPlane(1,3)] = max(FocusInfo(f).IntDenPerPlane(:,3));
else
end
if Channels>3
    [FocusInfo(f).MaxValue(1,4) FocusInfo(f).MaxPlane(1,4)] = max(FocusInfo(f).IntDenPerPlane(:,4));
else
end

end