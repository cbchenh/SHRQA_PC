%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial Traversing                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - Input parameters:                                                    %
%         img: img data (n * n * 1 or n * n * 3)                         %
%         l: level of the Hilbert Curve                                  %
%         v: orientation of the Hilbert Curve                            %
% - Output:                                                              %
%         ST: the Spatial Traversing Sequence of the input img           %
%                                                                        %
% Authors:                                                               %
%         1. Cheng-Bang Chen    email: cbchen@chengbangchen.me           %
%         2. Hui Yang           email: huy25@psu.edu                     %
%         3. Soundar Kumara     email: skumara@psu.edu                   %
%         Copyright 2019, Cheng-Bang Chen, All rights reserved.          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ST = SpatialTraverse(img, l, v)
    n1 = size(img, 1);
    n2 = size(img, 2);
    n3 = size(img, 3);
    if n1 ~= n2
       error('not a squared image.') 
    end
    img = int16(img); % Change image to int16 format
    img = imresize(img,[2^l 2^l]);
    M = Hilbertcurve(l,v,0);
    ST = zeros(0,n3);
    for i = 1:(2^(2*l))
        ST(i,:)=reshape(img(M(i,1),M(i,2),:),[1 n3]);
    end
end