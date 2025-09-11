set(0, 'DefaultFigureWindowStyle', 'docked')
%%
clear; close all;clc

PhotoDir = 'C:\Users\coxde\OneDrive\Masters\StaticTyreTestRig';

Files = dir(PhotoDir);
Files = Files(~[Files.isdir]);
FileNames = {Files.name};
FileNames = FileNames(1:end-3);


sigma = 5;
N = 100;

FinalImages = zeros(N, N, length(FileNames));

for i = 1:length(FileNames)
    fullName = fullfile(PhotoDir, FileNames{i});
    im = imread(fullName);
    
    % Convert to grayscale
    im_gray = rgb2gray(im);
    
    % Crop image to desired size
    im_cropped = imcrop(im_gray);

    imshow(im_cropped)

    % Apply gaussian blur before downsampling
    im_blurred = imgaussfilt(im_cropped, sigma);

    % Downsample to size NxN
    im_downsampled = imresize(im_blurred, [N, N]);

    imshow(im_downsampled)

    % Normalise values to range [0, 1]
    im_normalised = 1 - im2double(im_downsampled);

    % Save final image
    FinalImages(:, :, i) = im_normalised;
end

%%
figure
T = tiledlayout('flow');
T.TileSpacing = 'tight';
T.Padding = "compact";

for i = 1:length(FileNames)
    nexttile
    imshow(FinalImages(:, :, i))
end

%% 

X_lengths = [70, 85, 119.6, 134.7, 146.5, 166, 196, 198, 205];
Y_lengths = [68.9, 68.9, 104.9, 115.7, 122.3, 148, 164, 168, 175];
Z_loads = [150, 200, 250, 300, 350, 400, 500, 550, 600] * 9.81;

A_mat = [Z_loads(:), Z_loads(:).^2, log10(Z_loads(:))];

% Ax = b
X_vec_x = A_mat \ X_lengths(:);
X_vec_y = A_mat \ Y_lengths(:);

figure
hold on
scatter(Z_loads, X_lengths)
plot(Z_loads, A_mat * X_vec_x)
scatter(Z_loads, Y_lengths)
plot(Z_loads, A_mat * X_vec_y)
grid on
legend('Measured X', 'X Model', 'Measured Y', 'Y Model')
xlabel('Vertical Load [N]')
ylabel('Length [mm]')