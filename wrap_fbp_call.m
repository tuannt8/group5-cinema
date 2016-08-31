function rec = wrap_fbp_call( sinogram, crop_size, padsize_cor, padsize_roi, source_origin, detector_origin, angles, pixel_size, do_gpu, filter_type)

%% Crop sinogram
sinogram_cropped = crop_sinogram(sinogram, crop_size);

%% Pad sinogram for ROI data correction
sinogram_padded = pad_sinogram_roi(sinogram_cropped, padsize_roi);

%% Pad sinogram for COR correction
sinogram_ready = pad_sinogram_cor(sinogram_padded, padsize_cor);

%% Display ready sinogram
% figure
% show_image(sinogram_ready)
% colorbar
% colormap gray
% drawnow

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

vol_geom = astra_create_vol_geom(N,N);

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

%% Reconstruct using FBP
if do_gpu
    rec = fbp_gpu(-log(sinogram_ready), vol_geom, proj_geom, filter_type);
else
    rec = fbp_cpu(-log(sinogram_ready), proj_id);
end

end