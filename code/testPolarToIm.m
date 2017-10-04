%check polar to im 

support = 0:.01:1;

[A B] = meshgrid(support)

imSizeC = (length(support)*2);
imSizeD = (length(support)*2)-1;

C = PolarToIm(B,0,1,imSize,imSizeC);
D = PolarToIm(B,0,1,imSize,imSizeD);


figure;
subplot(1,3,1)
imagesc(B)
title('Grid B')
axis square
subplot(1,3,2)
imagesc(C)
title('PolarToIm with size')
axis square
subplot(1,3,3)
imagesc(D)
title('PolarToIm with size -1')
axis square