% Init
clear;
clc;

% Read image
image = imread("D:\EqualizeHistogram\EqualizeHistogram\image.png");

% Separate channels
redChannel = image(:,:,1);
greenChannel = image(:,:,2);
blueChannel = image(:,:,3);

% Equalized every histogram
equalizedRed = histeq(redChannel);
equalizedGreen = histeq(greenChannel);
equalizedBlue = histeq(blueChannel);

% Merge channels
equalizedImage = cat(3, equalizedRed, equalizedGreen, equalizedBlue);

% Plot
figure;

subplot(1, 2, 1);
imshow(image);
title('Original Image');

subplot(1, 2, 2);
imshow(equalizedImage);
title('Equalized Image');