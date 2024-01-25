%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularized Heterogeneous Recurrence Analysis Quantification           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Description:                                                         %
%      Input a sequence of integers, then provide a set of HRQAs &       %
%      corresponding IFS address.                                        %
% - Input:                                                               %
%         s: a sequence of integers                                      %
%         K: number of states (default: max(s))                          %
%         r: zoom in order (starts from 0, default: 1)                   %
%         a: iterated parameter (default:0.9999*sin(pi/K)/(1+sin(pi/K))) %
%         B: # of sections for calculating entropy & Gini (default: 20)  %
% - Output:                                                              %
%         IdxM: labels of states                                         %
%         HRR: Heterogeneous Recurrence Rate                             %
%         HMean: Heterogeneous Recurrence Mean                           %
%         HVar: Heterogeneous Recurrence Variance                        %
%         HSkew: Heterogeneous Recurrence Skewness                       %
%         HKurt: Heterogeneous Recurrence Kurtosis                       %
%         HEnt: Heterogeneous Recurrence Entropy                         %
%         HGini: Heterogeneous Recurrence Gini Index                     %
%                                                                        %
% Authors:                                                               %
%         Cheng-Bang Chen    email: cbchen@chengbangchen.me              %
%         Copyright 2021, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = RHRQA(s, K, r, a, B)
    if isempty(s)
        error('The input sequence, s, is empty.')
    end
    s = reshape(s,[length(s),1]);
    % Default values of K, r, a, B
    if nargin < 5
        if nargin < 4
           if nargin < 3
              if nargin < 2
                 K = max(s); 
              end
              r = 1;
           end
           a = 0.9999*sin(pi/K)/(1+sin(pi/K));
        end
        B=20;
    end
    if K < max(s)
       error('Error: K < max(s)');
    end
    % IFS address: Cv
    Cv = IFS(s, K, a);
    Cv = Cv(2:end,:);
    % Generate combinations of all local clusters: M
    IdxM = [];
    if(r == 0)
        IdxM = 0;
        Idx_set = 1:length(s);
        N = size(s,1);
        H_bar =[];
        for i = 1:K
            H_bar = [H_bar , sum(s==i)];
        end
        HRR = sum((H_bar./N).^2);
        LC_D = ((1/a)^(r))*pdist(Cv(Idx_set,:));
        HMean = mean(LC_D);
        HVar = (LC_D-HMean*ones(size(LC_D)))*(LC_D-HMean*ones(size(LC_D)))'/nchoosek(length(s),2); 
        HSkew = (sum((LC_D-HMean*ones(size(LC_D))).^3)/nchoosek(length(s),2))/(HVar)^(3/2); 
        HKurt = (sum((LC_D-HMean*ones(size(LC_D))).^4)/nchoosek(length(s),2))/(HVar)^(4/2);
        tmp_pb = [];
        tmp_pb = hist(LC_D, B);
        tmp_pb = nonzeros(tmp_pb);
        tmp_pb = tmp_pb./sum(tmp_pb);
        logpb = log2(tmp_pb);
        HEnt = -1*(tmp_pb'*logpb);
        HGini = 1 - (tmp_pb'*tmp_pb);
    else
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
        % Generate real local clusters
        LC = [];
        for i = 1:r
            LC(:,i)=s(i:end-r+i);
        end
        N = size(s,1); % The length of overall sequence
        % Unique local clusters
        UC = unique(LC, 'row');
        % Calculate the HRQAs for each unique local cluster
        % US_Res is the object to save the HRQAs for each local cluster:
        % [col1, col2, col3, col4, col5, col6, col7] =
        % [HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini]
        UC_Res = NaN(length(UC),7); 
        for i = 1:size(UC,1)
            %Idx_set = find(any(LC==UC(i,:),2));
            Idx_set = find(ismember(LC,UC(i,:),'row'));
            Idx_set = Idx_set + (r-1)*ones(size(Idx_set));
            % H_bar: Cardinality of local cluster
            H_bar = length(Idx_set);
            % HRR
            UC_Res(i,1) = (H_bar/N)^2;
            LC_ifs = Cv(Idx_set,:);
            LC_D = [];
    %         LC_D = pdist(LC_ifs);
            LC_D = ((1/a)^(r))*pdist(LC_ifs);
            if H_bar > 1        
                % HMean
                UC_Res(i,2) = mean(LC_D);
                % HVar
                UC_Res(i,3) = (LC_D-UC_Res(i,2)*ones(size(LC_D)))*(LC_D-UC_Res(i,2)*ones(size(LC_D)))'/nchoosek(H_bar,2);
                % HSkew
                UC_Res(i,4) = (sum((LC_D-UC_Res(i,2)*ones(size(LC_D))).^3)/nchoosek(H_bar,2))/(UC_Res(i,3))^(3/2);
                % HKurtosis
                UC_Res(i,5) = (sum((LC_D-UC_Res(i,2)*ones(size(LC_D))).^4)/nchoosek(H_bar,2))/(UC_Res(i,3))^(4/2);
                % HEnt
                % separate the distances of local cluster into B boxes & calculate
                % the probability of each box
                tmp_pb = [];
                tmp_pb = hist(LC_D, B);
                tmp_pb = nonzeros(tmp_pb);
                tmp_pb = tmp_pb./sum(tmp_pb);
                logpb = log2(tmp_pb);
                UC_Res(i,6) = -1*(tmp_pb'*logpb);

                % HGini
                UC_Res(i,7) = 1 - (tmp_pb'*tmp_pb);
            end
        end
        M_Res = NaN(size(IdxM,1),7);
        for i = 1:length(UC)
            ins_idx = find(ismember(IdxM,UC(i,:),'row'));
            if length(ins_idx)>0
                M_Res(ins_idx,:) = UC_Res(i,:);
            end
        end
        M_Res(isnan(M_Res(:,1)),1)=0;
        HRR = M_Res(:,1);
        HMean = M_Res(:,2); 
        HVar = M_Res(:,3); 
        HSkew = M_Res(:,4); 
        HKurt = M_Res(:,5); 
        HEnt = M_Res(:,6);
        HGini = M_Res(:,7);        
    end
end