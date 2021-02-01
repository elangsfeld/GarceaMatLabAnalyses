clc; clear all;

%DAPI = imread('C:\Users\Doug\Desktop\PDM Test\C1-Slide 09, RA, 0 Chase, FISH, 002_Reconstructed-1.tif');
TAg = imread('C:\Users\Doug\Desktop\PDM Test\C2-Slide 09, RA, 0 Chase, FISH, 002_Reconstructed-1.tif');
EdU = imread('C:\Users\Doug\Desktop\PDM Test\C4-Slide 09, RA, 0 Chase, FISH, 002_Reconstructed-1.tif');
Rows = size(TAg,1); Columns = size(TAg,2); Pixels = Rows*Columns;
PDM = zeros(Rows,Columns);
PDM_pos = zeros(Rows,Columns);
PDM_neg = zeros(Rows,Columns);
ZeroMask = zeros(Rows,Columns);

%DAPI2 = im2double(DAPI);
TAg2 = im2double(TAg);
EdU2 = im2double(EdU);

TAg_bg = imopen(TAg2,strel('disk',100));
EdU_bg = imopen(EdU2,strel('disk',100));


TAg3 = TAg2 - TAg_bg;
EdU3 = EdU2 - EdU_bg;

%bw = imbinarize(DAPI2); bw = bwareaopen(bw,10);
%bw = imerode(imdilate(imfill(bw,'holes'),strel('disk',50)),strel('disk',50));

%TAg_nuc = TAg3(bw); TAg_mean = mean(TAg_nuc);
%EdU_nuc = EdU3(bw); EdU_mean = mean(EdU_nuc);

TAg3_mean = reshape(TAg3,[Pixels,1]); TAg3_mean = mean(TAg3_mean);
EdU3_mean = reshape(EdU3,[Pixels,1]); EdU3_mean = mean(EdU3_mean);

for y = 1:Rows
       for x = 1:Columns
           TAg_DM(y,x) = TAg3(y,x) - TAg3_mean;
           EdU_DM(y,x) = EdU3(y,x) - EdU3_mean;
           PDM(y,x) = TAg_DM(y,x)*EdU_DM(y,x);
           
            if PDM(y,x) > 0
            PDM_pos(y,x) = PDM(y,x);
        else PDM_pos(y,x) == 0;
        end
        if PDM(y,x) < 0
            PDM_neg(y,x) = PDM(y,x);
        else PDM_neg(y,x) == 0;
        end
           
           %if TAg_DM(y,x)<0
            %   if EdU_DM(y,x)<0
             %      PDM(y,x)==0;
              % else
               %end
           %else PDM(y,x) = TAg_DM(y,x)*EdU_DM(y,x);
           %end
       end
end

for y = 1:Rows
    for x = 1:Columns
       if TAg_DM(y,x) < 0
           if EdU_DM(y,x) < 0
               ZeroMask(y,x) = 0;
           else ZeroMask(y,x) = 1;
           end
       else ZeroMask(y,x) = 1;
       end
    end
end

ZeroMask = logical(ZeroMask);
PDM_pos2 = PDM_pos.*ZeroMask;
PDM_neg2 = PDM_neg.*ZeroMask;

PDM_pos16 = im2uint16(PDM_pos2); imwrite(PDM_pos16, 'PDMpos.tiff','Compression','none');
PDM_neg16 = im2uint16(-PDM_neg2); imwrite(PDM_neg16, 'PDMneg.tiff','Compression','none');


%PDM = im2uint16(PDM);
%imwrite(PDM,'PDM.tiff','Compression','none');