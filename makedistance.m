function makedistance()
% Distance to vessel WALL (meters) from a fundus mask.
% Output: distance_to_vessel_m.csv  [x(m)  y(m)  ell(m)]
%
% Domain expected by COMSOL:
%   x in [0, W], y in [0, H], center at (W/2, H/2) = (0.03, 0.03)

close all; clc;

%% SETTINGS (to match COMSOL geometry)
fn = '27_manual1.gif';        % image (binary/near-binary); a dialog opens if not found
frame = 1;                    % first frame for GIFs

W = 565*0.000022;                     % physical width  [m]
H = 584*0.000022;                     % physical height [m]
% -> center is (0.03, 0.03)

% Image orientation fixes (apply BEFORE processing)
rotateSteps = 0;              % 0,1,2,3 => rotate by 0°,90°,180°,270° CCW
mirrorX_pre = false;          % mirror horizontally BEFORE processing
mirrorY_pre = false;          % mirror vertically   BEFORE processing

% Export orientation so CSV already matches COMSOL (origin bottom-left):
mirrorX_forCOMSOL = false;    % set true if you still see left-right inverted in COMSOL
mirrorY_forCOMSOL = true;     % usually true (fixes upside-down)

% Binarization/cleanup
minObjPixels   = 8;           % remove tiny specks
closingDiskRad = 1;           % closes micro-gaps (pixels)
smoothSigma_px = 0;           % 0=off; try 1..2 if you want very light smoothing

%% READ IMAGE
if ~isfile(fn)
    [f,p] = uigetfile({'*.gif;*.tif;*.tiff;*.png;*.jpg;*.jpeg;*.bmp','Images'}, ...
                      'Select vessel image'); 
    assert(~isequal(f,0),'No file selected.');  fn = fullfile(p,f);
end

[A,map] = imread(fn, frame);
if ~isempty(map)
    I = ind2gray(A,map);
else
    if ndims(A)==3, I = rgb2gray(A); else, I = A; end
    I = im2double(I);
end
I = mat2gray(I);

% Apply pre-orientation
if rotateSteps~=0, I = rot90(I, rotateSteps); end
if mirrorX_pre,     I = fliplr(I);           end
if mirrorY_pre,     I = flipud(I);           end

[ny,nx] = size(I);
px2m_x  = W / nx;         % meters per pixel (x)
px2m_y  = H / ny;         % meters per pixel (y)
fprintf('Scale: px2m_x=%.6g, px2m_y=%.6g m/pixel\n', px2m_x, px2m_y);

%% BINARIZE (Otsu) + auto-invert if necessary
th = graythresh(I);
BW = I > th;                       % assume bright vessels
if mean(BW(:)) > 0.5, BW = ~BW; end % if mostly white, invert (vessels were dark)

% light cleanup (NO hole fill, NO skeleton — we want WALL distance)
BW = bwareaopen(BW, minObjPixels);
BW = imclose(BW, strel('disk', closingDiskRad));

figure; imshow(BW); title('Vessel mask (white = vessel)');

%% EUCLIDEAN DISTANCE TO VESSEL WALL (in meters)
anis = abs(px2m_x - px2m_y)/max(px2m_x,px2m_y);
if anis < 0.02
    D_px = bwdist(BW, 'euclidean');
    D_m  = D_px * px2m_x;
else
    scale_y = px2m_x/px2m_y;
    ny_iso  = max(1, round(ny*scale_y));
    BW_iso  = imresize(BW, [ny_iso, nx], 'nearest');
    D_px_iso = bwdist(BW_iso,'euclidean');
    D_m_iso  = D_px_iso * px2m_x;
    D_m      = imresize(D_m_iso, [ny, nx], 'bilinear');
end

if smoothSigma_px > 0
    D_m = imgaussfilt(D_m, smoothSigma_px);
end

figure; imagesc(D_m); axis image off; colorbar; title('Distance to nearest vessel [m]');

%% BUILD COORDS IN METERS (origin bottom-left) + FINAL ORIENTATION FIX
% Create x,y axes (meters)
x_m = ((1:nx)-0.5) * px2m_x;        % 0..W
y_m = ((1:ny)-0.5) * px2m_y;        % top->bottom initially

% Flip to bottom->top so that (0,0) is bottom-left for COMSOL
D_out = D_m(end:-1:1,:);            % vertical flip
% Optional mirrors to match what you see in COMSOL
if mirrorY_forCOMSOL, D_out = D_out(end:-1:1,:); end  % extra vertical flip if needed
if mirrorX_forCOMSOL, D_out = fliplr(D_out);     end  % horizontal mirror if needed

% Recompute y-axis to be bottom-up
y_m = ((1:ny)-0.5) * px2m_y;        % 0..H bottom-up

[X_m, Y_m] = meshgrid(x_m, y_m);    % center is (W/2, H/2) = (0.03, 0.03)

%% WRITE CSV (ready for COMSOL Interpolation(File))
T = [X_m(:), Y_m(:), D_out(:)];
writematrix(T, 'distance_to_vessel_m.csv');
fprintf('Wrote: distance_to_vessel_m.csv (x[m], y[m], ell[m])\n');

%% QUICK OVERLAY CHECK
step_m = 10 * mean([px2m_x, px2m_y]);   % contour every ~10 px
levels = 0:step_m:max(D_out(:));
figure;
%imagesc(x_m, y_m, I(end:-1:1,:)); 
axis image; set(gca,'YDir','normal'); colormap(gray);
hold on; contour(x_m, y_m, D_out, levels, 'r', 'LineWidth', 0.8);
title('Overlay (origin bottom-left, center at 0.03,0.03)'); hold off;
end