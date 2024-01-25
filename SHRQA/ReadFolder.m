%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read all files from a folder                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Input parameters:                                                    %
%         folder_to_search: provide the directory                        %
%         filetype: provide the file type or pattern '*.png'             %
%         folder_to_output: provide the output directory                 %
%         l: the Hilbert Curve order                                     %
%         v: the Hilbert Curve orientation                               %
% Authors:                                                               %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         2. Hui Yang           email: huy25@psu.edu                     %
%         3. Soundar Kumara     email: skumara@psu.edu                   %
%         Copyright 2019, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
folder_to_search = 'C:\Users\CHENG-BANG\Dropbox\ImageData\Healthcare\BreastHistopathologyImages\prilimary\testing10';
filetype = fullfile(folder_to_search, '*.png');
folder_to_output = 'C:\Users\CHENG-BANG\Dropbox\ImageData\Healthcare\BreastHistopathologyImages\HRQA';
l = 6; % Hilbert Curve order
v = 1; % Hilbert Curve orientation
filelist = dir(filetype);
ST = [];
for Readinfile = 1 : length(filelist)
  File2Read = fullfile(folder_to_search, filelist(Readinfile).name);
  img_tmp = imread(File2Read); 
  ST_tmp = SpatialTraverse(img_tmp, l, v);
  ST = vertcat(ST,ST_tmp);
end
[Idx, Ub, Lb] = HAS(ST, size(ST,1)/16);
PlotCell(Ub, Lb);
for Readinfile = 1 : length(filelist)
  File2Read = fullfile(folder_to_search, filelist(Readinfile).name);
  img_tmp = imread(File2Read); 
  ST_tmp = SpatialTraverse(img_tmp, l, v);
  [s, ~] = SymbG(ST_tmp, Ub, Lb);
  [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = HRQA(s, size(Ub,1),1);
  Tmp_HRQA = horzcat(IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini);
  filename = filelist(Readinfile).name;filename = filename(1:end-4);
  File2Write = fullfile(folder_to_output, strcat(filename,'.csv'));
  csvwrite(File2Write,Tmp_HRQA);
end