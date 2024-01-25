%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Haar 2d wavelet                                                         %
%                                                                         %
% Before running the code:                                                %
% - set the working directory to the image folder                         %
%       [\your directory\all]                                             %
% - create 4 output folders to get wavelet results                        %
%       [\your directory\all_1A]                                          %
%       [\your directory\all_1D]                                          %
%       [\your directory\all_1H]                                          %
%       [\your directory\all_1V]                                          %
%                                                                         %
% - Cheng-Bang Chen: cxc1920@miami.edu                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
WD = 'D:\WD\all';
L = 1;
listing = dir(WD);

Dec={'A','H','V','D'};

for i = 1:L
    for j = 1:4
        folderName = strcat(WD,'_',num2str(i),Dec{j});
        if ~exist(folderName, 'dir')
            mkdir(folderName);
            disp('Folder created.');
        end
    end
end
tic
for i = 3:length(listing)
    name = strcat('D:\WD\all\',listing(i).name);
    sash = strfind(name,'\');
    outputname = name(sash(end)+1:end-4);
        img = imread(name);
        [A,H,V,D] = haart2(img,1,'integer');
        save(strcat('D:\WD\all_1A\',outputname,'.mat'),'A');
        save(strcat('D:\WD\all_1H\',outputname,'.mat'),'H');
        save(strcat('D:\WD\all_1V\',outputname,'.mat'),'V');
        save(strcat('D:\WD\all_1D\',outputname,'.mat'),'D');
    i
end
toc
