clear all; clc;
tic
cd(userpath);
%% Editables

Folder = 'FILEPATH/';
%Where are the images located? This filepath should end with a folder name, not a file name.
%IMPORTANT: Use format 'FILEPATH/'. The apostrophes and ending slash are important.

IgnoreChannels = 0; %Are there channels in the raw images that you want to ignore when creating image stacks? (1=yes, 0=no)
    IgChannels = [2]; %If yes (1), which channels do you want to ignore? If no (0), this input will be ignored.
%Ex: If you want to ignore channel 2, the input should be [2].
%Ex: If you want to ignore channels 2 and 4, the input should be [2 4];
%Ex: If you want to ignore channels 2, 4, and 5, the input should be [2 4 5];
%%%%%%%%%%%%%
%NOTE: You can use this script to save any number of channels you want, but
%       the HCS_ImageAnalysis script is only compatible with 2-4 channels.
%%%%%%%%%%%%%

%% Collection of Image Information %%
cd(Folder);
srcFiles = dir('*.tiff');
for u = 1:length(srcFiles)
    srcFiles_HiddenFilt(u,1) = srcFiles(u).name(1) == '.'; %Identifies the "hidden" files that Mac creates when transferring data between Windows and Macs.
end
srcFiles = srcFiles(~srcFiles_HiddenFilt); %Filters out these hidden files.

START = 1;
FINISH = length(srcFiles);

for f = START:FINISH
    if f == 1, disp('Collecting image information...'); else end;
    [pathstr,name,ext] = fileparts(srcFiles(f).name);
    
    location = strfind(name,'-');
    FileInfo(f,:).PlateLocation = name(1:12);
    FileInfo(f,:).Row = str2double(name(2:3));
    FileInfo(f,:).Column = str2double(name(5:6));
    FileInfo(f,:).Field = str2double(name(8:9));
    FileInfo(f,:).RCF = name(1:9);
    FileInfo(f,:).Plane = str2double(name(11:12));
    FileInfo(f,:).Channel = str2double(name(16));
    FileInfo(f,:).Remainder = name(17:end);
    Channel_Idx(f,1) = FileInfo(f).Channel;
end

%% Channel Indexing (only used if "Ignore Channels" = 1) %%

if IgnoreChannels ==1,
    for f = START:FINISH,
        for s = 1:numel(IgChannels),
            Ignore(f,s) = FileInfo(f).Channel == IgChannels(s);
        end
        Ignore(f,(numel(IgChannels)+1)) = sum(Ignore(f,:));
        Ignore(f,(numel(IgChannels)+2)) = ~Ignore(f,(numel(IgChannels)+1));
    end
    
    FileInfo = FileInfo(Ignore(:,(numel(IgChannels)+2)));
    Channel_Idx = Channel_Idx(Ignore(:,(numel(IgChannels)+2)));
    srcFiles = srcFiles(Ignore(:,(numel(IgChannels)+2)));
    FINISH = numel(FileInfo);
else
end

%% Image Collation and Saving %%

Channel_IdxUnique = unique(Channel_Idx);

for f = START:FINISH
    
    time(f,1).ElapsedSeconds = toc;
    clc
    filename = strcat(Folder,srcFiles(f).name);
    progress = ((FINISH-(FINISH-f))/FINISH)*100;
    progress2 = sprintf('On image %d of %d; %0.2f%c complete.',f,FINISH,progress,'%');
    disp(progress2);
    if progress < 10,
        disp('Estimated time remaining will display after 10% of images are analyzed...');
    else
        time(f).AverageSecondsPerLoop = time(f).ElapsedSeconds/(FINISH-(FINISH-f));
        time(f).EstimatedTotalSeconds = time(f).AverageSecondsPerLoop*(FINISH);
        time(f).EstimatedSecondsLeft = time(f).EstimatedTotalSeconds-time(f).ElapsedSeconds;
        time(f).EstimatedMinutesLeft = time(f).EstimatedSecondsLeft/60;
        time(f).EstimatedMinutesElapsed = time(f).ElapsedSeconds/60;
        estimate = sprintf('Run time: %0.2f minutes.',time(f).EstimatedMinutesElapsed);
        estimate2 = sprintf('Estimated time remaining: %0.2f minutes.',time(f).EstimatedMinutesLeft);
        disp(estimate);
        disp(estimate2);
    end
    cd(Folder);
    if f == 1
        mkdir('ImageStacks');          
    else
    end
        
    if f == 1 || size(strfind(FileInfo(f).RCF,FileInfo(f-1).RCF),1) > 0
        
        I = imread(srcFiles(f).name);
        for c = 1:numel(Channel_IdxUnique),
            if FileInfo(f).Channel == Channel_IdxUnique(c),
               IMAGE(:,:,FileInfo(f).Plane,c) = I; 
            else end;
        end
    
    else
        cd ./ImageStacks;
        disp('Saving image stack...');
        SAVE = strcat(FileInfo(f-1).RCF,'.ome.tiff');
        bfsave(IMAGE(:,:,:,:),SAVE);
        cd ../;
        I = imread(srcFiles(f).name);
        clearvars IMAGE;
        IMAGE(:,:,FileInfo(f).Plane,1) = I;
    end
    
    if f == FINISH
        cd ./ImageStacks;
        SAVE = strcat(FileInfo(f).RCF,'.ome.tiff');
        bfsave(IMAGE,SAVE);
        cd ../;
    else end    
end

toc