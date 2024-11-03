% Clear
clear;
clc;

% Step 1: Read image and convert to HSI
image = imread("D:\EqualizeHistogram\EqualizeHistogram\image.png");
[rows, cols, ~] = size(image);

% Convert RGB to HSI
tempImage = double(image) / 255;
hsiImage = rgb2hsi(tempImage);

% Extract H, S, and I channels
H = hsiImage(:, :, 1);
S = hsiImage(:, :, 2);
I = hsiImage(:, :, 3);

% Step 2: Set levels for intensity and saturation
L = 256;
M = 256;

% Initialize 2D histogram for intensity and saturation
histo2D = zeros(L, M);

% Populate the 2D histogram
for r = 1:rows
    for c = 1:cols
        iLevel = round(I(r, c) * (L - 1)) + 1;
        sLevel = round(S(r, c) * (M - 1)) + 1;
        histo2D(iLevel, sLevel) = histo2D(iLevel, sLevel) + 1;
    end
end

% Update parameters
n1_I = sum(sum(histo2D == 1, 2) > 0);
n2_I = sum(sum(histo2D == 2, 2) > 0);
b_I = n1_I / (n1_I + 2 * n2_I);

n1_S = sum(histo2D == 1, 'all');
n2_S = sum(histo2D == 2, 'all');
b_S = n1_S / (n1_S + 2 * n2_S);

% Step 3: Compute marginal CDF for intensity (F(x_I))
pdfS = zeros(M, 1);

for j = 1:M
    pdfS(j) = sum(histo2D(:, j)) / (rows * cols);
end

pdfI = zeros(L, 1);
cdfI = zeros(L, 1);

for i = 1:L
    if sum(histo2D(i, :) > 0)
        pdfI(i) = (sum(histo2D(i, :)) - b_I) / (rows * cols);
    else
        unseen_count = sum(sum(histo2D(:, :), 2) == 0);
        pdfI(i) = b_I * (L - unseen_count) / (rows * cols) / (unseen_count + eps);
    end
    
    if i > 1
        cdfI(i) = cdfI(i-1) + pdfI(i);
    else
        cdfI(i) = pdfI(i);
    end
end

pdfS_given_I = zeros(L, M);
cdfS_given_I = zeros(L, M);

for i = 1:L
    if cdfI(i) > 0
        sum_pdfS_unseen = sum(pdfS(histo2D(i, :) == 0));

        for j = 1:M
            if histo2D(i, j) > 0
                pdfS_given_I(i, j) = (histo2D(i, j) - b_S) / sum(histo2D(i, :));
            else
                unseen_count = sum(histo2D(i, :) == 0);
                pdfS_given_I(i, j) = b_S * (M - unseen_count) * pdfS(j) / sum(histo2D(i, :)) / sum_pdfS_unseen;
            end
            if j > 1
                cdfS_given_I(i, j) = cdfS_given_I(i, j - 1) + pdfS_given_I(i, j);
            else
                cdfS_given_I(i, j) = pdfS_given_I(i, j);
            end
        end
    end
end

% Step 4: Equalize I and S channels
equalizedI = zeros(rows, cols);
equalizedS = zeros(rows, cols);

for r = 1:rows
    for c = 1:cols
        % Original intensity and saturation levels
        iLevel = round(I(r, c) * (L - 1)) + 1;
        sLevel = round(S(r, c) * (M - 1)) + 1;
        
        % Compute the equalized values using the CDF
        equalizedI(r, c) = cdfI(iLevel);
        equalizedS(r, c) = cdfS_given_I(iLevel, sLevel);
    end
end

% Scale equalized intensity and saturation back to [0, 1] range
equalizedI = min(max(equalizedI, 0), 1);
equalizedS = min(max(equalizedS, 0), 1);

% Step 5: Reconstruct equalized HSI image and convert back to RGB
equalizedHSI = cat(3, H, equalizedS, equalizedI);
equalizedRGB = hsi2rgb(equalizedHSI);

equalizedRGB = uint8(equalizedRGB * 255);
equalizedRGB = double(equalizedRGB) / 255;

processedImage = zeros(rows, cols, 3);

% Calculate r, g, b
max_R = max(max(equalizedRGB(:, :, 1)));
max_G = max(max(equalizedRGB(:, :, 2)));
max_B = max(max(equalizedRGB(:, :, 3)));

for r = 1:rows
    for c = 1:cols
        % Get current RGB
        x_R = equalizedRGB(r, c, 1);
        x_G = equalizedRGB(r, c, 2);
        x_B = equalizedRGB(r, c, 3);
         
        r_ratio = x_R / max_R;
        g_ratio = x_G / max_G;
        b_ratio = x_B / max_B;
        l_x = r_ratio + g_ratio + b_ratio;
   
        % Calculate S-shape value
        s_shape_value = s_shape_transformation(l_x);
        
        % Judgement
        if s_shape_value / l_x <= 1
            processedImage(r, c, 1) = x_R * (s_shape_value / l_x);
            processedImage(r, c, 2) = x_G * (s_shape_value / l_x);
            processedImage(r, c, 3) = x_B * (s_shape_value / l_x);
        else
            % Convert RGB to CMY
            C = 1 - x_R;
            M = 1 - x_G;
            Y = 1 - x_B;
            
            factor = (3 - s_shape_value) / (3 - l_x);
            C_new = C * factor;
            M_new = M * factor;
            Y_new = Y * factor;
            
            % Convert CMY to RGB
            processedImage(r, c, 1) = (1 - C_new);
            processedImage(r, c, 2) = (1 - M_new);
            processedImage(r, c, 3) = (1 - Y_new);
        end
    end
end

% 3-D RGB (a)
R = zeros(rows * cols, 1);
G = zeros(rows * cols, 1);
B = zeros(rows * cols, 1);

index = 1;
for r = 1:rows
    for c = 1:cols
        R(index) = equalizedRGB(r, c, 1);
        G(index) = equalizedRGB(r, c, 2);
        B(index) = equalizedRGB(r, c, 3);
        index = index + 1;
    end
end

figure;
scatter3(R, G, B, 10, [R, G, B], 'filled');
xlabel('Red Channel');
ylabel('Green Channel');
zlabel('Blue Channel');
title('3D Scatter Plot of RGB Values');
grid on;
axis equal;

% 3-D RGB (b)
R = zeros(rows * cols, 1);
G = zeros(rows * cols, 1);
B = zeros(rows * cols, 1);

index = 1;
for r = 1:rows
    for c = 1:cols
        R(index) = processedImage(r, c, 1);
        G(index) = processedImage(r, c, 2);
        B(index) = processedImage(r, c, 3);
        index = index + 1;
    end
end

figure;
scatter3(R, G, B, 10, [R, G, B], 'filled');
xlabel('Red Channel');
ylabel('Green Channel');
zlabel('Blue Channel');
title('3D Scatter Plot of RGB Values');
grid on;
axis equal;

% Display results
figure;
subplot(1, 2, 1);
imshow(equalizedRGB);
title('Original Image');

subplot(1, 2, 2);
imshow(processedImage);
title('Equalized Image');

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

function y = s_shape_transformation(x)
    delta1 = 0;
    delta2 = 3;
    m = 1.5;
    n = 2;

    if x >= delta1 && x <= m
        y = delta1 + (m - delta1) * ((x - delta1) / (m - delta1))^n;
    elseif x > m && x <= delta2
        y = delta2 - (delta2 - m) * ((delta2 - x) / (delta2 - m))^n;
    else
        error('Input x is out of the defined range [delta1, delta2].');
    end
end