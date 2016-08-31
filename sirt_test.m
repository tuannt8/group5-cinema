%% Beware: Clear existing variable
clear
astra_clear

%% Parameters to specify
row_num = 512;

num_iter = 15; % Number of recon alg iterations to run

% Number of pixels to crop off each side of projections
crop_size = 0;

% Number of pixels to pad on each side projections for ROI correction
padsize_roi = 0;

% Number of pixels to pad asymmetrically for center-of-rotation correction,
% plus/minus specifies different sides.
padsize_cor = 4;

% Use GPU (true) or CPU (false)
do_gpu = true;

% Set filter for FBP - only works on GPU. Possible values:
% none, ram-lak, shepp-logan, cosine, hamming, hann, tukey, lanczos,
% triangular, gaussian, barlett-hann, blackman, nuttall, blackman-harris,
% blackman-nuttall, flat-top, kaiser, parzen
%filter_type = 'hamming';
filter_type = 'hamming';

%% Path to stored slices of reconstruction supplied by vendor

path_vendor_recon = ...
    '/work1/jakj/cinemax2/data/CaseGFRP/LFoV/Tomograms/Zeiss-bin/';
filename_vendor_str = ...
    'Butterfly_test_01_LFOV_50kV_VE_recon_Export%04d.bin';

angles_file = ...
    '/work1/jakj/cinemax2/data/CaseGFRP/LFoV/angles_LFoV.txt';
geometry_file = ...
    '/work1/jakj/cinemax2/data/CaseGFRP/LFoV/geometry_LFoV.txt';


%% Load sinogram and meta data and vendor reconstruction

% Load from MAT-file from ex03
load(sprintf('coursesinos/LFoV_sino%04d.mat',row_num));

% Read the angles from textfile
angles = importdata(angles_file);

% Extract the geometry parameters from textfile
[source_origin, detector_origin, pixel_size, proj_dims, ...
    proj_num, dims_vendor, recon_slice_num] = ...
    load_geometry(geometry_file);

% Load the corresponding vendor reconstruction
vendor_offset = (proj_dims(2)-recon_slice_num)/2;
vendor_row_num = row_num - vendor_offset;
fullfilename_vendor = fullfile(path_vendor_recon, ...
    sprintf(filename_vendor_str,vendor_row_num));
rec_vendor = load_recon_slice(fullfilename_vendor,dims_vendor);

%% Discard last projection and angle since copy of first
angles = angles(1:proj_num-1);
sinogram = sinogram(1:proj_num-1,:);

%% Display vendor reconstruction
figure
show_image(rec_vendor)
% axis([510 550 460 480])

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

% Number of pixels in object N-by-N
N = size(sinogram,2);

% Geometric magnification: Scaling factor mapping one object pixel in the
% center-of-rotation plane on to one detector pixel.
magnification = (source_origin+detector_origin)/source_origin;

% Size of a detector pixel relative to object pixel
detector_pixel_size = magnification;

% Number of detector pixels in ready sinogram
num_detector_pixels = size(sinogram_ready,2);

% Projections angles in radians (need to reverse by putting minus)
angles_astra = -angles/180*pi;

% Distance source to center-of-rotation, given in numbers of pixels.
source_center_dist = source_origin ./ pixel_size;

% % Distance detector to center-of-rotation, given in numbers of pixels.
detector_center_dist = detector_origin ./ pixel_size;

%% Set up ASTRA volume geometry

geom_padding = 32;

vol_geom = astra_create_vol_geom(N+geom_padding,N+geom_padding);

%% Set up ASTRA projection geometry

proj_geom = astra_create_proj_geom(...
    'fanflat', ...
    detector_pixel_size, ...
    num_detector_pixels, ...
    angles_astra,...
    source_center_dist,...
    detector_center_dist);

%% Set up ASTRA projector required for reconstruction on CPU

proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

%% Display geometry
% figure
% disp_geometry(proj_geom, vol_geom);
% title('Fan-beam scan geometry')

%% Load the sinogram data into a 2d data structure
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, -log(sinogram_ready));

%% Create a data object for the reconstruction
reconstruction_id = astra_mex_data2d('create', '-vol', vol_geom);

%% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);

%% Create configuration structure for reconstruction algorithm on GPU.
% Available algorithms: SIRT_CUDA, SART_CUDA, CGLS_CUDA, FBP_CUDA
cfg = astra_struct('CGLS_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;

%% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

%% Run num_iter iterations of the algorithm

num_iters = (10:10:150);

solving_time = zeros(1, numel(num_iters));

% reconstructions = zeros(numel(num_iters), N, N);

% reconstructions = cell(numel(num_iters,1));

% for i = numel(num_iters)
    tic
%     astra_mex_algorithm('iterate', alg_id, num_iters(5));
    astra_mex_algorithm('iterate', alg_id, num_iter);
%     solving_time(i) = toc;
    solving_time = toc;
%     reconstructions(i, :, :) = astra_mex_data2d('get', rec_id);
    reconstruction = astra_mex_data2d('get', rec_id);
%     reconstructions{i} = astra_mex_data2d('get', rec_id);
%     break;
% end

% tic
% astra_mex_algorithm('iterate', alg_id, num_iter);
% solving_time = toc;

% disp(['Solution time for running SIRT:' num2str(solving_time)]);

%% Get the result
% rec = astra_mex_data2d('get', rec_id);

%% Display

% figure
% imagesc(rot90(rec,2));
% axis image off
% colorbar
% title_str = sprintf('N = %d',N);
% title(title_str)

% for i = numel(num_iters)
    figure
%     imshow(reconstruction, []);
    show_image(reconstruction);
%     imagesc(reconstruction);
%     axis image off
    colorbar;
    title('SIRT, num iter = ');
% end

    tic
%     astra_mex_algorithm('iterate', alg_id, num_iters(5));
    astra_mex_algorithm('iterate', alg_id, num_iter);
%     solving_time(i) = toc;
    solving_time2 = toc;

    reconstruction2 = astra_mex_data2d('get', rec_id);

        figure
%     imshow(reconstruction, []);
    show_image(reconstruction2);
%     imagesc(reconstruction);
%     axis image off
    colorbar;
    title('SIRT, num iter = ');
