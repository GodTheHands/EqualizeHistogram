% Init
clear;
clc;

% Read image
imageRGB = imread("D:\EqualizeHistogram\EqualizeHistogram\image.png");
[rows, cols, ~] = size(imageRGB);

% Convert RGB to HSI
imageRGB = double(imageRGB) / 255; % Normalize RGB values to [0,1]
imageHSI = rgb2hsi(imageRGB);

% Extract H, S, and I channels
H = imageHSI(:, :, 1);
S = imageHSI(:, :, 2);
I = imageHSI(:, :, 3);

% Define the PDF for the intensity channel
L = 256;
intensityValues = linspace(0, 1, L);

pdfI = zeros(1, L);
for k = 1:L
    if intensityValues(k) <= 0.5
        pdfI(k) = 12 * intensityValues(k)^2;
    else
        pdfI(k) = 12 * (1 - intensityValues(k))^2;
    end
end
pdfI = pdfI / L * 2; % Normalize PDF

% Calculate the CDF from the PDF
cdfI = cumsum(pdfI);

% Perform histogram equalization on the intensity channel
I_eq = zeros(rows, cols);
for r = 1:rows
    for c = 1:cols
        % Find the nearest index in intensityValues to the pixel's intensity
        [~, idx] = min(abs(intensityValues - I(r, c)));
        
        % Map intensity value using the computed CDF
        I_eq(r, c) = cdfI(idx);
    end
end

% Construct the equalized HSI image
equalizedHSI = cat(3, H, S, I_eq);

% Convert equalized HSI back to RGB
equalizedImageRGB = hsi2rgb(equalizedHSI);
equalizedImageRGB = uint8(equalizedImageRGB * 255);

% Plot the original and equalized images
figure;
subplot(1, 2, 1);
imshow(imageRGB);
title('Original Image in RGB');

subplot(1, 2, 2);
imshow(equalizedImageRGB);
title('Equalized Image in RGB');

% Supporting function: RGB to HSI
function hsi = rgb2hsi(rgb)
    R = rgb(:, :, 1);
    G = rgb(:, :, 2);
    B = rgb(:, :, 3);
    numerator = 0.5 * ((R - G) + (R - B));
    denominator = sqrt((R - G).^2 + (R - B).*(G - B));
    theta = acos(numerator ./ (denominator + eps));
    H = theta;
    H(B > G) = 2 * pi - H(B > G);
    H = H / (2 * pi);
    S = 1 - 3 * min(cat(3, R, G, B), [], 3) ./ (R + G + B + eps);
    I = (R + G + B) / 3;
    hsi = cat(3, H, S, I);
end

% Supporting function: HSI to RGB
function rgb = hsi2rgb(hsi)
    H = hsi(:, :, 1) * 2 * pi;
    S = hsi(:, :, 2);
    I = hsi(:, :, 3);
    R = zeros(size(H));
    G = zeros(size(H));
    B = zeros(size(H));
    idx = (0 <= H) & (H < 2 * pi / 3);
    B(idx) = I(idx) .* (1 - S(idx));
    R(idx) = I(idx) .* (1 + S(idx) .* cos(H(idx)) ./ cos(pi / 3 - H(idx)));
    G(idx) = 3 * I(idx) - (R(idx) + B(idx));
    idx = (2 * pi / 3 <= H) & (H < 4 * pi / 3);
    H(idx) = H(idx) - 2 * pi / 3;
    R(idx) = I(idx) .* (1 - S(idx));
    G(idx) = I(idx) .* (1 + S(idx) .* cos(H(idx)) ./ cos(pi / 3 - H(idx)));
    B(idx) = 3 * I(idx) - (R(idx) + G(idx));
    idx = (4 * pi / 3 <= H) & (H <= 2 * pi);
    H(idx) = H(idx) - 4 * pi / 3;
    G(idx) = I(idx) .* (1 - S(idx));
    B(idx) = I(idx) .* (1 + S(idx) .* cos(H(idx)) ./ cos(pi / 3 - H(idx)));
    R(idx) = 3 * I(idx) - (G(idx) + B(idx));
    rgb = cat(3, R, G, B);
end
