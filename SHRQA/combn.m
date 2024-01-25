%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get all combinations of the sequences                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Description:                                                         %
%      Input the cardinality and the length of sequence to get all       %
%      combinations.                                                     %
% - Input:                                                               %
%         K: the number of cardinality                                   %
%         r: length of sequence                                          %
% - Output:                                                              %
%         IdxM: a list of all combinations                               %
%                                                                        %
% Authors:                                                               %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         2. Hui Yang           email: huy25@psu.edu                     %
%         3. Soundar Kumara     email: skumara@psu.edu                   %
%         Copyright 2019, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function IdxM = combn(K,r)
    IdxM = [];
    for i = 1:r
        tmp_c = [];
        tmp_e = [];
        for j = 1:K
            tmp_e = [tmp_e;j*ones(K^(r-i),1)];
        end
        for k = 1:(K^(i-1))
            tmp_c = [tmp_c;tmp_e]; 
        end
        IdxM(:,r-i+1)=tmp_c;
    end
end