%% Beware: Clear existing variable
clear
astra_clear

%% Parameters to specify
row_num = 512;

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

projection_path = ...
    '/work1/jakj/cinemax2/data/CaseGFRP/LFoV/Projections/bin/';
projection_file_str = ...
    'Butterfly_test_01_LFOV_50kV_VE_Export%04d.bin';

%% Read projections and create 3d sinogram
slice_num = 5;
slice_range = (-slice_num:1:slice_num) + row_num;
num_projections = 801;
detector_size = [1024,1024];

sinogram = zeros(801,1024,numel(slice_range));

for i = 1:num_projections
    disp(['Loading projection i = ' num2str(i)]);
    projection = load_projection([projection_path sprintf(projection_file_str,i)], detector_size);
    % Size 801 x 1024 x #slice_range:
    sinogram(i, :, :) = projection(slice_range, :)';
end

%% Load sinogram and meta data and vendor reconstruction

% Load from MAT-file from ex03
% load(sprintf('coursesinos/LFoV_sino%04d.mat',row_num));

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
% Remove the last projection (overlapping with the first)
sinogram = sinogram(1:proj_num-1,:,:);

%% Display vendor reconstruction
% figure
% show_image(rec_vendor)
% axis([510 550 460 480])

%% Crop sinogram
% 3d cropping?
% sinogram_cropped = crop_sinogram(sinogram, crop_size);

%% Pad sinogram for ROI data correction
% 3d padding?
% sinogram_padded = pad_sinogram_roi(sinogram_cropped, padsize_roi);

%% Pad sinogram for COR correction
% sinogram_ready = pad_sinogram_cor(sinogram_padded, padsize_cor);

sinogram_ready = pad_sinogram_cor_3d(sinogram, padsize_cor);

%% Display ready sinogram (of center slice)
figure
% show_image(sinogram_ready(:,:,11))
show_image(sinogram(:,:,6))
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
num_detector_pixels = size(sinogram,2);

% Projections angles in radians (need to reverse by putting minus)
angles_astra = -angles/180*pi;

% Distance source to center-of-rotation, given in numbers of pixels.
source_center_dist = source_origin ./ pixel_size;

% % Distance detector to center-of-rotation, given in numbers of pixels.
detector_center_dist = detector_origin ./ pixel_size;

%% Set up ASTRA volume geometry

% Start by just having a 800 x 1024 x #slice_range
vol_geom = astra_create_vol_geom(N,N,numel(slice_range));

% Create the data structure for reconstruction volume
rec_id = astra_mex_data3d('create', '-vol', vol_geom);

%% Set up ASTRA projection geometry

detector_row_count = numel(slice_range);
detector_column_count = num_detector_pixels;

proj_geom = astra_create_proj_geom( ...
    'cone', ...
    detector_pixel_size, ...
    detector_pixel_size, ...
    detector_row_count, ...
    detector_column_count, ...
    angles_astra, ...
    source_center_dist, ...
    detector_center_dist);

%% Load the sinogram data into a 3d data structure

% Rearrange data for the to fit with projection geometry
sinogram_id = astra_mex_data3d('create', '-sino', proj_geom, -log(permute(sinogram,[2,1,3])));

%% Display geometry
% figure
% disp_geometry(proj_geom, vol_geom);
% title('Fan-beam scan geometry')

%% Create configuration structure for reconstruction algorithm on GPU.
% Available algorithms: SIRT_CUDA, SART_CUDA, CGLS_CUDA, FBP_CUDA
cfg = astra_struct('FDK_CUDA');
% cfg = astra_struct('SIRT3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ShortScan = true;

%% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

%% Run num_iter iterations of the algorithm
disp('Performing reconstrution using FDK running for iterations')
tic
astra_mex_algorithm('iterate', alg_id);
run_time = toc();
disp(['Took ' num2str(run_time) ' to run algorithm'])

%% Retrieve our reconstruction
reconstruction = astra_mex_data3d('get', rec_id);
figure
for i = 1:size(reconstruction,3)
    show_image(rot90(squeeze(reconstruction(:,:,i))));
    drawnow
    pause(1)
end

%% Shows the FDK center slice for comparisonv

% This should show the standard center slice for comparison
figure
show_image(rot90(squeeze(reconstruction(:,:,6))));
