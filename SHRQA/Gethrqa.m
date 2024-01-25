addpath('C:\Users\cxc1920\Dropbox\matlab\SFC\SHRQA')

[rr,cc] = size(M);
hrqa = [];
for idx = 1:cc
    tic
    s = M(:,idx);
    [IdxM, HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini] = HRQA(s, length(C));
    Tmp_HRQA = vertcat(HRR, HMean, HVar, HSkew, HKurt, HEnt, HGini)';
    % hrqa(idx,:) = Tmp_HRQA;
    fname = sprintf('hrqa_%d.csv',idx)
    csvwrite(fname,Tmp_HRQA) 
    toc
end
% save HRQA.mat hrqa