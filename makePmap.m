%% makePmap+macula.m 
clear; close all; clc;

%% ===================== PARAMETRI (X MODIFICA) =====================
imgFile    = '27_manual1.gif';   % file immagine (vasi bianchi su nero)

% Dominio fisico (m)
W = 565*0.000022;  H = 584*0.000022;

% Orientamento immagine
rotateSteps = 0;      % 0/1/2/3 => 0°/90°/180°/270° CCW
mirrorX_pre = false;  % flip orizzontale
mirrorY_pre = false;  % flip verticale

% Migliorie per rami sottili / transizioni
up         = 2;       % upsample (1=off, 2 consigliato)
sig_skel   = 1.2;     % smoothing lungo scheletro [px]
sig_vessel = 0.8;     % smoothing guidato dentro vaso [px]

% Correzione raggio (PSF/segmentazione)
gammaR = 1.0;         % 0.7–1.2 se diametri paiono gonfi/magri

% Legge P ~ R^{-alpha}
P0   = 7e-8;          % [m/s]
Rref = 20e-6;         % [m]
alpha= 0.6;           % 0.3–0.5 fisiol, 0.8–1.0 marcato
Rmin = 4e-6;          % [m]
Pmin = 1e-8;          % [m/s]
Pmax = 1e-5;          % [m/s]

% Hotspot maculare (ellisse)
xc   = 0.03;  yc = 0.03;     % [m]
a    = 1.2e-2; b  = 0.9e-2;  % [m]
theta= 0;      sCNV = 0.10e-3;  % [rad], [m]
P_CNV= 3e-6;                  % [m/s]


%% Lettura immagine
if ~isfile(imgFile)
    [f,p] = uigetfile({'*.gif;*.tif;*.tiff;*.png;*.jpg;*.jpeg;*.bmp','Images'}, ...
                      'Seleziona immagine dei vasi');
    if isequal(f,0), error('Nessun file selezionato.'); end
    imgFile = fullfile(p,f);
end

[A,map] = imread(imgFile,1);
if ~isempty(map)
    I = ind2gray(A,map);
else
    if ndims(A)==3, I = rgb2gray(A); else, I = A; end
    I = im2double(I);
end
I = mat2gray(I);

% Orientamento pre
if rotateSteps~=0, I = rot90(I, rotateSteps); end
if mirrorX_pre,     I = fliplr(I); end
if mirrorY_pre,     I = flipud(I); end

% Upsample
if up>1, I = imresize(I, up, 'bicubic'); end

[ny,nx] = size(I);
px2m_x = W/nx;  px2m_y = H/ny;
fprintf('Scala: %.6g (x), %.6g (y) m/px — up=%d\n', px2m_x, px2m_y, up);

%% Segmentazione leggera (immagine già pulita)
th = graythresh(I);
BW = I > th; 
if mean(BW(:))>0.5, BW = ~BW; end
BW = bwmorph(BW,'bridge',Inf);
BW = imclose(BW, strel('disk',1));
figure; imshow(BW); title('Maschera vasi (bianco = vaso)');

%% RAGGIO LOCALE R(x,y) [m]  — versione robusta senza branching
% Rendiamo i pixel “isotropi” in y e calcoliamo sempre con lo stesso percorso
scale_y = px2m_x/px2m_y;
ny_iso  = max(1, round(ny*scale_y));
BW_iso  = imresize(BW, [ny_iso, nx], 'nearest');
D_in_m  = bwdist(~BW_iso,'euclidean') * px2m_x;        % distanza interna in m
D_in_m  = imresize(D_in_m, [ny, nx], 'bilinear');      % torna alla griglia originale

% Skeleton (senza togliere rami corti)
if exist('bwskel','file')
    Sk = bwskel(BW,'MinBranchLength',0);
else
    Sk = bwmorph(BW,'skel',Inf);
end

% Raggio su centerline + smoothing lungo scheletro
R_skel = zeros(size(D_in_m));  R_skel(Sk) = D_in_m(Sk);
if exist('imgaussfilt','file')
    num = imgaussfilt(R_skel, sig_skel);
    den = imgaussfilt(double(Sk), sig_skel) + eps;
else
    h = fspecial('gaussian', max(5,ceil(6*sig_skel)), sig_skel);
    num = imfilter(R_skel, h, 'replicate');
    den = imfilter(double(Sk), h, 'replicate') + eps;
end
R_skel = (num ./ den) .* double(Sk);

% Propagazione + smoothing guidato dentro vaso
[~,idx] = bwdist(Sk);
Rmap = gammaR * R_skel(idx);     % raggio ovunque [m]
M = double(BW);
if exist('imgaussfilt','file')
    num = imgaussfilt(Rmap.*M, sig_vessel);
    den = imgaussfilt(M,      sig_vessel) + eps;
else
    h = fspecial('gaussian', max(5,ceil(6*sig_vessel)), sig_vessel);
    num = imfilter(Rmap.*M, h, 'replicate');
    den = imfilter(M,      h, 'replicate') + eps;
end
Rmap = num ./ den;

%% Griglia COMSOL (origine in basso-sx)
x_axis = ((1:nx)-0.5)*px2m_x;    % 0..W
y_axis = ((1:ny)-0.5)*px2m_y;    % 0..H
R_out  = Rmap(end:-1:1,:);       % flip verticale
BW_out = BW(end:-1:1,:);
[Xm,Ym] = meshgrid(x_axis, y_axis);

%% P(R) + hotspot maculare
Rclip = max(R_out, Rmin);
Pbase = P0 .* (Rref ./ Rclip) .^ alpha;

xr =  (Xm - xc)*cos(theta) + (Ym - yc)*sin(theta);
yr = -(Xm - xc)*sin(theta) + (Ym - yc)*cos(theta);
chi = (xr./a).^2 + (yr./b).^2;
sigma = sCNV / max(a,b) + eps;
wCNV  = 1 ./ (1 + exp((chi - 1)/sigma));   % ~1 dentro ellisse

Pmap = Pbase.*(1 - wCNV) + P_CNV*wCNV;
Pmap = min(Pmax, max(Pmin, Pmap));

%% Salvataggi CSV
writematrix([Xm(:), Ym(:), R_out(:)], 'radius_map_m.csv');
writematrix([Xm(:), Ym(:), Pmap(:)],   'Pmap_mps.csv');
fprintf('Scritto:  radius_map_m.csv   e   Pmap_mps.csv\n');

%% Verifica #1: R solo sui bordi
B = bwperim(BW); 
B_out = B(end:-1:1,:);
R_border = R_out; R_border(~B_out) = NaN;
figure; imagesc(x_axis, y_axis, R_border);
axis image; set(gca,'YDir','normal'); colorbar;
title('R sul bordo del vaso [m] (smussato)');

% Overlay su fundus
If = I(end:-1:1,:);
figure; imagesc(x_axis, y_axis, If); colormap(gray); axis image; set(gca,'YDir','normal'); hold on;
h = imagesc(x_axis, y_axis, R_border); 
set(h, 'AlphaData', ~isnan(R_border)*0.85); colorbar;
title('Overlay: fundus + raggio sul bordo'); hold off;

%% Verifica #2: istogramma diametri
R_bordo    = R_out(B_out);
D_bordo_um = 2*R_bordo*1e6;
figure; histogram(D_bordo_um, 40);
xlabel('Diameter [\mum]'); ylabel('Count');
title('Distribution diameter on edges');
med = median(D_bordo_um,'omitnan'); 
q25 = prctile(D_bordo_um,25); 
q75 = prctile(D_bordo_um,75);
fprintf('Diametro mediana = %.1f µm (IQR %.1f–%.1f µm)\n', med, q25, q75);

%% Vasi colorati per P (CLASSI) — con maschera
Pv = Pmap(BW_out); Pv = Pv(:);
edges = prctile(Pv(~isnan(Pv)), [0 25 50 75 100]);  % quartili
Pcls = nan(size(Pmap));
for k=1:4
    maskK = BW_out & Pmap>=edges(k) & Pmap<=edges(k+1);
    Pcls(maskK) = k;
end

figure; 
imagesc(x_axis, y_axis, Pcls, 'AlphaData', double(BW_out)*0.95);
axis image; set(gca,'YDir','normal');
colormap([0 0 1; 0 0.8 1; 1 0.85 0; 1 0 0]);
colorbar('Ticks',1:4,'TickLabels', ...
   {sprintf('%.2g–%.2g',edges(1),edges(2)), sprintf('%.2g–%.2g',edges(2),edges(3)), ...
    sprintf('%.2g–%.2g',edges(3),edges(4)), sprintf('%.2g–%.2g',edges(4),edges(5))});
title('Coloured vessels for P (classes, masked)');
