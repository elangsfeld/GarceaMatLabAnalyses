

function [Ch1,Ch2,Ch3,Ch4,Focus_Ch1,Focus_Ch2,Focus_Ch3,Focus_Ch4,ZFocus] = InFocusImage(I,ZPlanes,ResY,ResX,Channels,f)

for m = 1:ZPlanes
    Ch1(:,:,m) = I{1,1}{m,1}(:,:);

    Ch2(:,:,m) = I{1,1}{ZPlanes+m,1}(:,:);
    
    if Channels>2, Ch3(:,:,m) = I{1,1}{(2*ZPlanes)+m,1}(:,:); 
    else Ch3 = zeros(1); end
    
    if Channels>3, Ch4(:,:,m) = I{1,1}{(3*ZPlanes)+m,1}(:,:); 
    else Ch4 = zeros(1); end
end

for m = 1:ZPlanes
    if sum(Ch1(:,:,m),'all') == 0 && m>1
        Ch1(:,:,m) = Ch1(:,:,(m-1));
    elseif sum(Ch1(:,:,m),'all') == 0 && m==1
        Ch1(:,:,m) = Ch1(:,:,(m+1));
    else   
        Ch1_grad(m).GradImage(:,:) = imgradient(Ch1(:,:,m));
        Ch1_grad(m).Mean = mean(Ch1_grad(m).GradImage(:,:),'all');
        Ch1_grad(m).Percentile = prctile(Ch1_grad(m).GradImage(:,:),[95 99],'all');
        Filtered(:,:,m) = Ch1_grad(m).GradImage(Ch1_grad(m).Percentile(1,1) < Ch1_grad(m).GradImage < Ch1_grad(m).Percentile(2,1));    
        FilteredMean(m,1) = mean(mean(Filtered(:,:,m)));
    end
    
    if sum(Ch2(:,:,m),'all') == 0 && m>1
        Ch2(:,:,m) = Ch2(:,:,(m-1));
    elseif sum(Ch2(:,:,m),'all') == 0 && m==1
        Ch2(:,:,m) = Ch2(:,:,(m+1));
    else end
    
    if Channels > 2
        if sum(Ch3(:,:,m),'all') == 0 && m>1
            Ch3(:,:,m) = Ch3(:,:,(m-1));
        elseif sum(Ch3(:,:,m),'all') == 0 && m==1
            Ch3(:,:,m) = Ch3(:,:,(m+1));
        else end
    else
    end
    
    if Channels > 3
        if sum(Ch4(:,:,m),'all') == 0 && m>1
            Ch4(:,:,m) = Ch4(:,:,(m-1));
        elseif sum(Ch4(:,:,m),'all') == 0 && m==1
            Ch4(:,:,m) = Ch4(:,:,(m+1));
        else end
    else
    end
end
    
    
    [MAX(f,1),ZFocus(f,1)] = max(FilteredMean);

if ZFocus(f,1) == 1 && ZPlanes > 3
    for y = 1:ResY
        for x = 1:ResX
            Focus_Ch1(y,x) = max([Ch1(y,x,ZFocus(f,1)),Ch1(y,x,(ZFocus(f,1)+1))]);
            Focus_Ch2(y,x) = max([Ch2(y,x,ZFocus(f,1)),Ch2(y,x,(ZFocus(f,1)+1))]);
            if Channels>2, Focus_Ch3(y,x) = max([Ch3(y,x,ZFocus(f,1)),Ch3(y,x,(ZFocus(f,1)+1))]); else end;
            if Channels>3, Focus_Ch4(y,x) = max([Ch4(y,x,ZFocus(f,1)),Ch4(y,x,(ZFocus(f,1)+1))]); else end;
        end
    end
elseif ZFocus(f,1) == ZPlanes && ZPlanes > 3
    for y = 1:ResY
        for x = 1:ResX
            Focus_Ch1(y,x) = max([Ch1(y,x,ZFocus(f,1)),Ch1(y,x,(ZFocus(f,1)-1))]);
            Focus_Ch2(y,x) = max([Ch2(y,x,ZFocus(f,1)),Ch2(y,x,(ZFocus(f,1)-1))]);
            if Channels>2, Focus_Ch3(y,x) = max([Ch3(y,x,ZFocus(f,1)),Ch3(y,x,(ZFocus(f,1)-1))]); else end;
            if Channels>3, Focus_Ch4(y,x) = max([Ch4(y,x,ZFocus(f,1)),Ch4(y,x,(ZFocus(f,1)-1))]); else end;
        end
    end
elseif ZPlanes > 3
    for y = 1:ResY
        for x = 1:ResX
            Focus_Ch1(y,x) = max([Ch1(y,x,ZFocus(f,1)),Ch1(y,x,(ZFocus(f,1)+1)),Ch1(y,x,(ZFocus(f,1)-1))]);
            Focus_Ch2(y,x) = max([Ch2(y,x,ZFocus(f,1)),Ch2(y,x,(ZFocus(f,1)+1)),Ch2(y,x,(ZFocus(f,1)-1))]);
            if Channels>2, Focus_Ch3(y,x) = max([Ch3(y,x,ZFocus(f,1)),Ch3(y,x,(ZFocus(f,1)+1)),Ch3(y,x,(ZFocus(f,1)-1))]); else end;
            if Channels>3, Focus_Ch4(y,x) = max([Ch4(y,x,ZFocus(f,1)),Ch4(y,x,(ZFocus(f,1)+1)),Ch4(y,x,(ZFocus(f,1)-1))]); else end;
        end
    end
elseif ZPlanes < 4
    Focus_Ch1 = max(Ch1,[],3);
    Focus_Ch2 = max(Ch2,[],3);
    if Channels>2, Focus_Ch3 = max(Ch3,[],3); else end
    if Channels>3, Focus_Ch4 = max(Ch4,[],3); else end
end

if exist('Focus_Ch3') == 0, Focus_Ch3 = Ch3; else end
if exist('Focus_Ch4') == 0, Focus_Ch4 = Ch4; else end

end