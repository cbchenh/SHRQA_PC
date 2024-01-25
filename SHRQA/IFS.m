%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterated Function System                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Description:                                                         %
%      Input a sequence of integers, then provide a set of IFS address   %
%      accordingly.                                                      %
% - Input:                                                               %
%         s: a sequence of integers                                      %
%         K: number of states (default: max(s))                          %
%         a: iterated parameter (default:0.9999*sin(pi/K)/(1+sin(pi/K))) %
% - Output:                                                              %
%         Cv: IFS addresses of s                                         %
%                                                                        %
% Authors:                                                               %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         2. Hui Yang           email: huy25@psu.edu                     %
%         3. Soundar Kumara     email: skumara@psu.edu                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Cv] = IFS(s, K, a, opt)
    if isempty(s)
        error('The input sequence, s, is empty.')
    end
    s = reshape(s,[length(s) 1]);
    if nargin < 4
        opt = 0;
    end
    if nargin < 3
        if nargin < 2
            K = max(s);
        end
        a = 0.9999*sin(pi/K)/(1+sin(pi/K));
    end
    if isempty(a)
       a = 0.9999*sin(pi/K)/(1+sin(pi/K)); 
    end
    Cv = [0, 0];
    for i = 1:length(s)
        Cv(i+1,:) = [a*Cv(i,1)+cos(2*pi*s(i)/K),a*Cv(i,2)+sin(2*pi*s(i)/K)];
    end
    if opt > 0
        figure();
        plot(Cv(:,1),Cv(:,2),'.');
        xlim([-2 2]);
        ylim([-2 2]);
        title('IFS address')
    end
end