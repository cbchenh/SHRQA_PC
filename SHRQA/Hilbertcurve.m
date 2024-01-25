%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hilbert Space-Filling Curve                                            %
% - Input:                                                               %
%         l: the level of Hilbert Curve                                  %
%         v: the orientation of the Hilbert Curve                        %
%            v={1,2,3,4} == {N, E, S, W}                                 %
%         opt: 0: no image; 1: plot image                                %
% - Output:                                                              %
%         The route of the Hilbert Curve desired                         %
%                                                                        %
% Author:                                                                %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         2. Hui Yang           email: huy25@psu.edu                     %
%         3. Soundar Kumara     email: skumara@psu.edu                   %
%         Copyright 2019, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = Hilbertcurve(l,v,opt)
    if nargin <3
        opt = 1;
    end
    % Initialization:
    N = zeros(0,2);
    E = zeros(0,2);
    S = zeros(0,2);
    W = zeros(0,2);
    % Generate 4 orientations of Hilbert Curves recurrsively
    for i = 1:l
        [N,E,S,W] = HC(N,E,S,W);
    end
    
    switch v
        case 1
            M = N;
        case 2
            M = E;
        case 3
            M = S;
        case 4
            M = W; 
        otherwise
            disp("Invalid input v!!")
            M = [];
    end
    
    M = [[0,0];M];
    M(:,1) = cumsum(M(:,1));
    M(:,2) = cumsum(M(:,2));
    minC1 = min(M(:,1));
    minC2 = min(M(:,2));
    for i = 1:size(M,1)
        M(i,1) = M(i,1) - minC1 +1;
        M(i,2) = M(i,2) - minC2 +1;
    end
    if opt ~= 0
        plot(M(:,1),M(:,2),'Linewidth',3)
        axis off
        axis equal
    end
end
%%
function [newN,newE,newS,newW] = HC(N,E,S,W)
    newN = [E;[0,1];N;[1,0];N;[0,-1];W];
    newE = [N;[1,0];E;[0,1];E;[-1,0];S];
    newS = [W;[0,-1];S;[-1,0];S;[0,1];E];
    newW = [S;[-1,0];W;[0,-1];W;[1,0];N];
end
