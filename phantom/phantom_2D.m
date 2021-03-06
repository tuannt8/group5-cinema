%% Create phantom data
clear;
clc;
close all;

vol_geom = astra_create_vol_geom(128, 128, 128);

angles = linspace2(0, 2*pi, 360);

%proj_geom = astra_create_proj_geom('parallel3d', 1.0, 1.0, 128, 192, angles);


% proj_geom = astra_create_proj_geom('cone', ...
%                                 1  % det_spacing_x, 
%                                 1  % det_spacing_y, 
%                                 det_row_count, ...
%                                 det_col_count, ...
%                                 angles, ...
%                                 source_origin, ...
%                                 origin_det)
proj_geom = astra_create_proj_geom('cone', ...
                                    1, ...
                                    1, ...
                                    128, ...
                                    128, ...
                                    angles, ...
                                    500, ...
                                    128);


% Create a simple hollow cube phantom
cube = zeros(128,128,128);
cube(17:112,17:112,17:112) = 1;
cube(33:96,33:96,33:96) = 0;

% Create projection data from this
[proj_id, proj_data] = astra_create_sino3d_cuda(cube, proj_geom, vol_geom);

% Add noise?

% Display a single projection image
figure, imshow(squeeze(proj_data(:,20,:))',[])

%% Extract one slide
% row_num = 50;
% psize = size(proj_data);
% sinogram = zeros(psize(1), psize(2));
% for i = 1:psize(2)
%     sinogram(:,i) = proj_data(:,i,row_num);
% end
% sinogram = sinogram';

row_num = 64;
sinogram = proj_data(:,:,row_num)';

%%
%% Parameters to specify
crop_size = 0;

padsize_roi = 0;

padsize_cor = 0;

% Use GPU (true) or CPU (false)
do_gpu = true;

% Set filter for FBP - only works on GPU
filter_type = 'shepp-logan';

%%
%% Crop sinogram
sinogram_cropped = crop_sinogram(sinogram, crop_size);

%% Pad sinogram for ROI data correction
sinogram_padded = pad_sinogram_roi(sinogram_cropped, padsize_roi);

%% Pad sinogram for COR correction
sinogram_ready = pad_sinogram_cor(sinogram_padded, padsize_cor);

%% Display ready sinogram
figure
show_image(sinogram_ready)
colorbar
colormap gray
drawnow

%% Extract parameters for ASTRA
source_origin = proj_geom.DistanceOriginSource;
detector_origin = proj_geom.DistanceOriginDetector;


%

% Number of pixels in object N-by-N
N = size(sinogram,2);

% Geometric magnification: Scaling factor mapping one object pixel in the
% center-of-rotation plane on to one detector pixel.
detector_pixel_size = (source_origin+detector_origin)/source_origin;

% Number of detector pixels in ready sinogram
num_detector_pixels = size(sinogram_ready,2);


%% Set up ASTRA volume geometry
%N = 128;
vol_geom = astra_create_vol_geom(N,N);

%% Set up ASTRA projection geometry

proj_geom = astra_create_proj_geom(...
    'fanflat', ...
    detector_pixel_size, ...
    num_detector_pixels, ...
    angles,...
    source_origin,...
    detector_origin);

%% Set up ASTRA projector required for reconstruction on CPU

proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

% Display geometry
figure
disp_geometry(proj_geom, vol_geom);
title('Fan-beam scan geometry')

% Reconstruct using FBP
if do_gpu
    rec = fbp_gpu(sinogram_ready, vol_geom, proj_geom, filter_type);
else
    rec = fbp_cpu(sinogram_ready, proj_id);
end

% FBP-reconstructed image
figure
show_image(rec)
title('FBP reconstruction')

