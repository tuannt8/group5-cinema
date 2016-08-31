%% Beware: Clear existing variables
clear
astra_clear

%% Parameters to experiment with for region-of-interest reconstruction

% Object size N-by-N pixels:
N = 512; 

% Number of detector pixels
num_detector_pixels = 128;

% Number of pixels to pad on each side of sinogram to fix ROI problem
padsize = 64;

% Use GPU (true) or CPU (false)
do_gpu = true;

% Set filter for FBP - only works on GPU
filter_type = 'Ram-Lak';

%% Fixed parameters

% Number of projection angles.
num_angles = 1024;

% Distance source to center-of-rotation, given in numbers of pixels.
source_origin = 1500;

% % Distance detector to center-of-rotation, given in numbers of pixels.
detector_origin = 500;

% Incident X-ray flux in photons per time unit. Lower=noisier.
flux = 1e6; 

%% Derived parameters

% Geometric magnification: Scaling factor mapping one object pixel in the
% center-of-rotation plane on to one detector pixel.
magnification = (source_origin+detector_origin)/source_origin;

% Size of a detector pixel relative to object pixel
detector_pixel_size = magnification;  

% Angles at which projections are taken. Here equally spaced, 360 degrees.
angles = linspace2(0,2*pi,num_angles);

%% Set up ASTRA volume geometry

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

%% Display geometry
fig_geom = figure;
clf
disp_geometry(proj_geom, vol_geom);
title('Fan-beam scan geometry')

%% Generate test image and sinogram data
X0 = phantom(N);
[~,sinogram] = astra_create_sino(X0, proj_id);

% Add noise
if ~isinf(flux)
    sinogram_noisy = astra_add_noise_to_sino(sinogram,flux);
else 
    sinogram_noisy = sinogram;
end

%% Display clean and noisy sinogram
figure
show_image(sinogram)
title('Clean sinogram data')

figure
show_image(sinogram_noisy)
title('Noisy sinogram data')

%% Reconstruct using FBP
if do_gpu
    rec = fbp_gpu(sinogram_noisy, vol_geom, proj_geom, filter_type);
else
    rec = fbp_cpu(sinogram_noisy, proj_id);
end

%% Display original and FBP-reconstructed images
figure
show_image(X0)
title('Original image')
caxis([0,1])
hold on
plot_roi(N,num_detector_pixels)
plot(N/2*[1,1],ylim,'-b')

figure
show_image(rec)
title('FBP reconstruction')
caxis([0,1])
hold on
plot_roi(N,num_detector_pixels)
plot(N/2*[1,1],ylim,'-b')

%% Until here nothing has changed from ex01. Now: sinogram padding

% Copy projection geometry and increase number of detectors for padding
proj_geom_pad = proj_geom;
proj_geom_pad.DetectorCount = num_detector_pixels + 2*padsize;

% Also update the projector with the new projection geometry
proj_id_pad = astra_create_projector('line_fanflat',proj_geom_pad,vol_geom);

% Pad the sinogram with copies of the left/right-most pixels
sinogram_noisy_pad = pad_sinogram_roi(sinogram_noisy,padsize);

% Do FBP reconstruction from padded data.
if do_gpu
    rec_pad = fbp_gpu(sinogram_noisy_pad, vol_geom, proj_geom_pad, filter_type);
else
    rec_pad = fbp_cpu(sinogram_noisy_pad, proj_id_pad);
end

%% Display

% Show padded geometry on top of unpadded
figure(fig_geom)
disp_geometry(proj_geom_pad, vol_geom);
title('Padded fan-beam scan geometry')

% Display the padded sinogram
figure
show_image(sinogram_noisy_pad)
title('Padded noisy sinogram')

% Display the FBP reconstruction from padded data
figure
show_image(rec_pad)
title('FBP reconstruction from padded data')
caxis([0,1])
hold on
plot_roi(N,num_detector_pixels)
plot(N/2*[1,1],ylim,'-b')

% Display line profiles of the original, FBP from unpadded and padded data
figure
clf
plot(X0(:,N/2))
hold on
plot(rec(:,N/2))
plot(rec_pad(:,N/2),'--')
plot(0.5+0.5*(N+num_detector_pixels*[-1,-1,nan,1,1]),[ylim,nan,ylim],'-k')
legend('Original', 'FBP', 'FBP padded','ROI boundaries')
title('Central vertical profiles')
