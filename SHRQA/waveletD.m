%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  wavelet decomposition for image
%  - Cheng-Bang Chen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H1img,V1img,D1img,A2img,H2img,V2img,D2img] = waveletD(img)
    [c,s] = wavedec2(img, 2, 'haar');
    [H1,V1,D1] = detcoef2('all',c,s,1);
    [H2,V2,D2] = detcoef2('all',c,s,2);
    A2 = appcoef2(c,s,'haar',2);
    V1img = wcodemat(V1,255,'mat',1);
    H1img = wcodemat(H1,255,'mat',1);
    D1img = wcodemat(D1,255,'mat',1);
    V2img = wcodemat(V2,255,'mat',2);
    H2img = wcodemat(H2,255,'mat',2);
    D2img = wcodemat(D2,255,'mat',2);
    A2img = wcodemat(A2,255,'mat',2);
end