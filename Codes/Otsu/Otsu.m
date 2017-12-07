I = imread('100007.jpg');
%I = imnoise(I,'salt & pepper');
tic
k = 1;
level = multithresh(rgb2gray(I), k);
subplot(1,3,1); imshow(I)
title('Original Image');
subplot(1,3,2); histogram(I, 'BinWidth', 0.001)
subplot(1,3,2); histogram(I, 'BinWidth', 1)
title('Histogram of Pixel Brightness');
xlabel('Pixel Brightness');
ylabel('Number of Pixels');
if(k == 1)
hold on; subplot(1,3,2); plot(double(200)*ones(1,2501), 0:1:2500);
end
otsu_L = imquantize(rgb2gray(I),level);

subplot(1,3,3); imagesc(otsu_L);
toc