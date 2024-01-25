%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symblic Sequence Generator                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Input:                                                               %
%         X: multichannel signals                                        %
%         Ub: upper bound list of all regions                            %
%         Lb: lower bound list of all regions                            %
% - Output:                                                              %
%         S: symbolic sequence                                           %
%         DistM: distance Matrix                                         %
% Authors:                                                               %
%         Cheng-Bang Chen    email: cbchen@chengbangchen.me              %
%         Hui Yang           email: huy25@psu.edu                        %
%         Soundar Kumara     email: skumara@psu.edu                      %
%         09/27/2019                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S, DistM] = SymbG(X, Ub, Lb)
    if nargin < 4
        if nargin < 3
            error('Please provide boundries of each region.')
        end
        if size(Ub) ~= size(Lb)
            error('The sizes of ub and lb are different.')
        end
        Cb = zeros(size(Ub));
        for i = 1:size(Ub,1)
            Cb(i,:) = mean([Ub(i,:);Lb(i,:)]);
        end
    end
    
    if size(Ub) ~= size(Cb)
        error('The sizes of ub and cb are different.')
    end
    
    S = zeros(size(X,1),1);
    for idx = 1:size(Ub,1)
        dummy_idx = 1:size(X,1);
        for d_id = 1:size(Ub,2)
            dummy_idx = intersect(dummy_idx, find(X(:,d_id)<Ub(idx,d_id) & X(:,d_id)>=Lb(idx,d_id)));
        end
        S(dummy_idx) = idx;
    end
    DistM = squareform(pdist(Cb));
end