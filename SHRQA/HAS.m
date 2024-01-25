%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hyperoctree Aggregate Segmentation                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Input:                                                               %
%         X: multichannel signals                                        %
%         Cap: capacity of one discretized state                         %
%-----------------------   Hidden Inputs   ------------------------------%
%         ub: upper bounds of the current region                         %
%         lb: lower bounds of the current region                         %
%         deepth: the deepth of search                                   %
%-----------------------   Hidden Inputs   ------------------------------%
% - Output:                                                              %
%         Idx: labels of segments                                        %
%         Ub: Upper bounds of segments                                   %
%         Lb: Lower bounds of segments                                   %
%                                                                        %
% Authors:                                                               %
%         Cheng-Bang Chen    email: cbchen@chengbangchen.me              %
%         Copyright 2019, Cheng-Bang Chen, All rights reserved.          %
% -- update (08/17/2021): deepth = 5 for default value                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Idx, Ub, Lb] = HAS(X, Cap, ub, lb, deepth)
    set(0, 'RecursionLimit', 100)    
    if nargin < 5
        deepth = 7;
        if nargin == 3 | nargin ==4
            error('Please check the upper bounds and lower bounds.')
        end
        if nargin < 3
            if nargin < 2
                error('Please check the input vector X and set 1Capacity.')
            else
                ub = max(X)+ 0.000001*(max(X)-min(X));
                lb = min(X)- 0.000001*(max(X)-min(X));
            end
        end
    end
    
    cb = mean([ub;lb]);
    Ub = [];
    Lb = [];
    % n: # of objects in the current region
    % d: # of dimensions of the region
    [n, d]  = size(X);
    % if n > Cap, then segment to 2^d subspaces.
    if (n > Cap) & (deepth > 0)
        newub = vertcat(ub,cb);
        newlb = vertcat(cb,lb);
        IdxM = combn(2,d);
        [iter, chk] = size(IdxM);
        for k = 1:iter
            nub = zeros(1,chk);
            nlb = zeros(1,chk);
            dummy_idx = 1:size(X,1);
            for d_id = 1:chk
                nub(d_id) = newub(IdxM(k,d_id),d_id);
                nlb(d_id) = newlb(IdxM(k,d_id),d_id);
                dummy_idx = intersect(dummy_idx, find(X(:,d_id)<nub(d_id) & X(:,d_id)>=nlb(d_id)));
            end
            tmpX = X(dummy_idx,:);
            [~, tmpUb, tmpLb] = HAS(tmpX, Cap, nub, nlb, (deepth-1));
            Ub = vertcat(Ub,tmpUb);
            Lb = vertcat(Lb,tmpLb);
            tmpUb =[];
            tmpLb =[];            
        end
    else
        Ub = ub;
        Lb = lb;
    end
    Idx = 1:size(Ub,1);
end