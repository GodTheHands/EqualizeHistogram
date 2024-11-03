% Init
clear;
clc;

% Read image
image = imread("D:\EqualizeHistogram\EqualizeHistogram\image.png");
filteredImage = imgaussfilt(image, 1.5);
image = filteredImage;

[rows, cols, ~] = size(image);

% Set gray levels
L = 256;

% Calculate 3-D histogram
histo3D = zeros(L, L, L);
for r = 1:rows
    for c = 1:cols
        R = image(r, c, 1) + 1;
        G = image(r, c, 2) + 1;
        B = image(r, c, 3) + 1;
        histo3D(R, G, B) = histo3D(R, G, B) + 1;
    end
end

histo3D = histo3D ./ (rows * cols);

% Calculate joint CDF
cdf3D = cumsum(cumsum(cumsum(histo3D, 1), 2), 3);

% Calculate uniformly distributed CDF
uniformCDF = zeros(L, L, L);
for k = 0:L-1
    for s = 0:L-1
        for t = 0:L-1
            uniformCDF(k + 1, s + 1, t + 1) = (k + 1) * (s + 1) * (t + 1) / L^3;
        end
    end
end

% Initialize equalized image
equalizedImage = zeros(rows, cols, 3, 'uint8');
% Choose smallest tuple
for r = 1:rows
    for c = 1:cols
        R = image(r, c, 1) + 1;
        G = image(r, c, 2) + 1;
        B = image(r, c, 3) + 1;
        
        % Get original CDF
        y = cdf3D(R, G, B);
        
        % Find tuple s.t. constraints
        kPrime = R - 1;
        sPrime = G - 1;
        tPrime = B - 1;
        
        found = false;
        lastKPrime = kPrime;
        lastSPrime = sPrime;
        lastTPrime = tPrime;
        
        while y ~= 0
            currentUniformCDF = uniformCDF(kPrime + 1, sPrime + 1, tPrime + 1);
            
            while currentUniformCDF >= y
                lastKPrime = kPrime;
                lastSPrime = sPrime;
                lastTPrime = tPrime;
                
                if kPrime > 0
                    kPrime = kPrime - 1;
                end
                
                if sPrime > 0
                    sPrime = sPrime - 1;
                end
                
                if tPrime > 0
                    tPrime = tPrime - 1;
                end

                currentUniformCDF = uniformCDF(kPrime + 1, sPrime + 1, tPrime + 1);
                found = true;
            end
            
            if found == true
                equalizedImage(r, c, 1) = lastKPrime;
                equalizedImage(r, c, 2) = lastSPrime;
                equalizedImage(r, c, 3) = lastTPrime;
                break;
            end
            
            while currentUniformCDF < y && ~(kPrime == L - 1 && sPrime == L - 1 && tPrime == L - 1)
                if kPrime < L - 1
                    kPrime = kPrime + 1;
                end
                
                if sPrime < L - 1
                    sPrime = sPrime + 1;
                end
                
                if tPrime < L - 1
                    tPrime = tPrime + 1;
                end
                
                currentUniformCDF = uniformCDF(kPrime + 1, sPrime + 1, tPrime + 1);
            end
            
            equalizedImage(r, c, 1) = kPrime; 
            equalizedImage(r, c, 2) = sPrime;
            equalizedImage(r, c, 3) = tPrime;
            break;
        end
    end
end

% Plot
figure;
subplot(1, 2, 1);
imshow(image);
title('Original Image');

subplot(1, 2, 2);
imshow(equalizedImage);
title('Equalized Image');