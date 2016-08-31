%% Beware: Clear existing variables
clear
astra_clear

%% Parameters to specify

% Object size N-by-N pixels:
N = 512; 

% Number of detector pixels
num_detector_pixels = 512;

% Number of projection angles.
num_angles = 1024;

% Distance source to center-of-rotation, given in numbers of pixels.
source_origin = 1500;

% % Distance detector to center-of-rotation, given in numbers of pixels.
detector_origin = 500;

% Incident X-ray flux in photons per time unit. Lower=noisier.
flux = 1e6; 

% Use GPU (true) or CPU (false)
do_gpu = true;

% Set filter for FBP - only works on GPU. Possible values:
% none, ram-lak, shepp-logan, cosine, hamming, hann, tukey, lanczos,
% triangular, gaussian, barlett-hann, blackman, nuttall, blackman-harris,
% blackman-nuttall, flat-top, kaiser, parzen
filter_type = 'shepp-logan';

%% Derived parameters

% Geometric magnification: Scaling factor mapping one object pixel in the
% center-of-rotation plane on to one detector pixel.
magnification = (source_origin+detector_origin)/source_origin;

% Size of a detector pixel relative to object pixel
detector_pixel_size = magnification;  

% Angles at which projections are taken. Here equally spaced, 360 degrees.
angles = linspace2(0,2*pi,num_angles);

%% Set up ASTRA volume geometry

% vol_geom = astra_create_vol_geom(row_count, col_count);
%
%  Create a 2D volume geometry.  See the API for more information.
%  row_count: number of rows.
%  col_count: number of columns.
%
%  vol_geom: MATLAB struct containing all information of the geometry.

vol_geom = astra_create_vol_geom(N,N);

%% Set up ASTRA projection geometry

% proj_geom = astra_create_proj_geom('fanflat', det_width, det_count, ...
%                       angles, source_origin, origin_det)
% 
%  Create a 2D flat fan beam geometry.  See the API for more information.
% 
%  det_width: distance between two adjacent detectors
%  det_count: number of detectors in a single projection
%  angles: projection angles in radians, should be between -pi/4 and 7pi/4
%  source_origin: distance between the source and the center of rotation
%  origin_det: distance between the center of rotation and the detector array
%
%  proj_geom: MATLAB struct containing all information of the geometry

proj_geom = astra_create_proj_geom(...
    'fanflat', ...
    detector_pixel_size, ...
    num_detector_pixels, ...
    angles,...
    source_origin,...
    detector_origin);

%% Set up ASTRA projector required for reconstruction on CPU

%  proj_id = astra_create_projector(type, proj_geom, vol_geom)
%  
%  Create a new projector object based on projection and volume geometry.  
%  Used when the default values of each projector are sufficient.  
% 
%  type: type of the projector.  'line_fanflat', 'strip_fanflat'   See API for more information.
%  proj_geom: MATLAB struct containing the projection geometry.
%  vol_geom: MATLAB struct containing the volume geometry.
%
%  proj_id: identifier of the projector as it is now stored in the astra-library.
  
proj_id = astra_create_projector('line_fanflat', proj_geom, vol_geom);

%% Display geometry
figure
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

t = 0;

%% Reconstruct using FBP
if do_gpu
    tic
    rec = fbp_gpu(sinogram_noisy, vol_geom, proj_geom, filter_type);
    t = toc();
    disp('Running on GPU');
else
    tic
    rec = fbp_cpu(sinogram_noisy, proj_id);
    t = toc();
    disp('Running on CPU') ;
end

disp(['Time taken to run: ' num2str(t)]);

%% Display original and FBP-reconstructed images
figure
show_image(X0)
title('Original image')

figure
show_image(rec)
title('FBP reconstruction')

figure
show_image(abs(X0-rec));
title('Difference image')
colorbar;