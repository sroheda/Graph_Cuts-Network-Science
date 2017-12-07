tic
I = imread('100007.jpg');
I = imnoise(I, 'salt & pepper');
A = rgb2gray(I);
A = im2double(A); %A./max(max(max(A)));
G=fspecial('gaussian', [50 50], 5); 
Ig=imfilter(A,G,'same'); 
JK=A-Ig; imshow(JK);
toc