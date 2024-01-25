%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHRQA                                                                   %
%                                                                         %
% Before running the code:                                                %
% - create a folder to store the results                                  %
%       [\your directory\output]                                          %
%                                                                         %
% - Cheng-Bang Chen: cxc1920@miami.edu                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
clc
WD = 'D:\WD\all';
L = 1;
addpath 'F:\Dropbox\matlab\SFC\SHRQA'
filelist = dir(WD);

tic
l = 8; % Hilbert Curve order
v = 1; % Hilbert Curve orientation
ST = [];
for Readinfile = 3:length(filelist)
  img_tmp = imread(strcat(filelist(Readinfile).folder,'\',filelist(Readinfile).name));
  if size(size(img_tmp),2)==2
      tmp = uint8(256*ones(2^l))-img_tmp;
      img_tmp = cat(3, tmp, tmp, tmp);
  end
  ST_tmp = SpatialTraverse(img_tmp, l, v);
  ST = vertcat(ST,ST_tmp);
  Readinfile
end
toc
tic
[Idx, Ub, Lb] = HAS(ST, size(ST,1)/8);
S = SymbG(ST, Ub, Lb);
toc

tic
output = [];
namelist = strings(length(filelist)-2,1);
for Readinfile = 3:length(filelist)
  File2Read = strcat(filelist(Readinfile).folder,'\',filelist(Readinfile).name);
  strt = strfind(File2Read,'\');
  namelist(Readinfile-2) = File2Read(strt(end)+1:end-4);
  img_tmp = imread(File2Read);
  if size(size(img_tmp),2)==2
      tmp = uint8(256*ones(2^l))-img_tmp;
      img_tmp = cat(3, tmp, tmp, tmp);
  end
  ST_tmp = SpatialTraverse(img_tmp, l, v);
  [s, ~] = SymbG(ST_tmp, Ub, Lb);
  %   IFS(s, size(Ub,1), 0.045 ,1);
  [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = RHRQA(s, size(Ub,1),1);
  Tmp_HRQA = [HRR', HMean', HVar', HSkew', HKurt', HEnt', HGini'];
  output(Readinfile-2,:) = Tmp_HRQA;
end
varname = ["HRR";"HMean";"HVar";"HSkew";"HKurt";"HEnt";"HGini"];
index = repelem(varname,length(HRR));
varname = index + "_" + repmat(1:length(HRR),1,7)';
output = array2table(output,'RowNames',namelist,'VariableNames',varname');
outname = strcat('D:\WD\output\gs_orig.csv');
writetable(output,outname,'WriteRowNames',true)
toc

clear
clc

addpath 'F:\Dropbox\matlab\SFC\SHRQA'
filelist = dir('D:\WD\');

l = 7; % Hilbert Curve order
v = 1; % Hilbert Curve orientation

parfor folder = 4:7 %% change 4:7 to your setting
    output = [];
    foldername = filelist(folder).name;
    listing = dir(strcat(filelist(folder).folder,'\',foldername));
    namelist = strings(length(listing)-2,1);
    tic
    ST=[];
    for Readinfile = 3:length(listing)
      File2Read = strcat(listing(Readinfile).folder,'\',listing(Readinfile).name);
      strt = strfind(File2Read,'\');
      namelist(Readinfile-2) = File2Read(strt(end)+1:end-4);
      img_tmp = load(File2Read);
      img_tmp = cell2mat(struct2cell(img_tmp));
      ST_tmp = SpatialTraverse(img_tmp, l, v);
      ST = vertcat(ST,ST_tmp);
    end
    toc

    tic
    [Idx, Ub, Lb] = HAS(ST, size(ST,1)/8);
    S = SymbG(ST, Ub, Lb);
    toc
    
    %save(strcat('F:\ProstateCancer\output\HAS',num2str(folder),'.mat'),"Idx","Ub","Lb","S");
    parsave(strcat('D:\WD\output\HAS',num2str(folder),'.mat'), Idx, Ub, Lb, S);
    tic
    for Readinfile = 3:length(listing)
      File2Read = strcat(listing(Readinfile).folder,'\',listing(Readinfile).name);
      strt = strfind(File2Read,'\');
      namelist(Readinfile-2) = File2Read(strt(end)+1:end-4);
      img_tmp = load(File2Read);
      img_tmp = cell2mat(struct2cell(img_tmp));
      ST_tmp = SpatialTraverse(img_tmp, l, v);
      [s, ~] = SymbG(ST_tmp, Ub, Lb);
      [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = RHRQA(s, size(Ub,1),1);
      Tmp_HRQA = [HRR', HMean', HVar', HSkew', HKurt', HEnt', HGini'];
      output(Readinfile-2,:) = Tmp_HRQA;
    end
    varname = ["HRR";"HMean";"HVar";"HSkew";"HKurt";"HEnt";"HGini"];
    index = repelem(varname,length(HRR));
    varname = index + "_" + repmat(1:length(HRR),1,7)';
    output = array2table(output,'RowNames',namelist,'VariableNames',varname');
    outname = strcat('D:\WD\output\',foldername((end-1):end),'.csv');
    writetable(output,outname,'WriteRowNames',true)
    toc
    fprintf(strcat("-----finish ",foldername((end-1):end),"-----\n"))
end
