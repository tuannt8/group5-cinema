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
padsize_cor = 2;

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

pad_size = (-10:1:10);


for i = 1:numel(pad_size)
    rec = wrap_fbp_call( sinogram, crop_size, pad_size(i), padsize_roi, source_origin, detector_origin, angles, pixel_size, do_gpu, filter_type);
    figure;
    show_image(rec);
    title(['FBP with pad size =' num2str(pad_size(i))]);
    pause(2);
end

%% FBP-reconstructed image
figure
show_image(rec)
title('FBP reconstruction')

shiftx = +9;
axis([510+shiftx 550+shiftx 460 480])
