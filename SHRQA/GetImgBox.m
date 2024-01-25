%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discretize State-space by giving sample data                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Input parameters:                                                    %
%         folderName: directory of Sampled Images                        %
%         imgType: provide the file type or pattern '*.png'              %
%         l: the Hilbert Curve order                                     %
%         v: the Hilbert Curve orientation                               %
%         BoxCap: the capacity for each subspace                         %
% - Output:                                                              %
%         Idx: the indexes of the discretized subspaces                  %
%         Ub: the Upperbound of each discretized subspace                %
%         Lb: the Lowerbound of each discretized subspace                %
% Authors:                                                               %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         2. Hui Yang           email: huy25@psu.edu                     %
%         3. Soundar Kumara     email: skumara@psu.edu                   %
%         Copyright 2019, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% folderName = 'C:\Users\CHENG-BANG\Dropbox\ImageData\Healthcare\BreastHistopathologyImages\prilimary\testing10';
% imgType = '*.png';
% l = 6; % Hilbert Curve order
% v = 1; % Hilbert Curve orientation
% BoxCap = 320000;
function [Idx, Ub, Lb] = GetImgBox(folderName, imgType, l, v, BoxCap)
    folder_to_search = folderName;
    filetype = fullfile(folder_to_search, imgType);
    filelist = dir(filetype);
    ST = [];
    for Readinfile = 1 : length(filelist)
      File2Read = fullfile(folder_to_search, filelist(Readinfile).name);
      img_tmp = imread(File2Read); 
      % img_tmp = imresize(img_tmp, [2^l 2^l]);
      ST_tmp = SpatialTraverse(img_tmp, l, v);
      ST = vertcat(ST, ST_tmp);
    end
    [Idx, Ub, Lb] = HAS(ST, BoxCap);
end