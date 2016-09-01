%% Create phantom data
clear;
clc;
close all;

% Create a simple hollow cube phantom
load phantom.mat;
cube = phantom(1:128,1:128,1:128);
s = size(cube);

vol_geom = astra_create_vol_geom(s(1), s(2), s(3));

angles = linspace2(0, 2*pi, 200);

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
                                    250, ...
                                    250, ...
                                    angles, ...
                                    300, ...
                                    60);




% Create projection data from this
[proj_id, proj_data] = astra_create_sino3d_cuda(cube, proj_geom, vol_geom);

% Add noise
sigma = 10;
proj_data = proj_data + sigma*randn(size(proj_data));

% Display a single projection image
figure, imshow(squeeze(proj_data(:,20,:))',[])

% Create a data object for the reconstruction
rec_id = astra_mex_data3d('create', '-vol', vol_geom);

% Set up the parameters for a reconstruction algorithm using the GPU
% SIRT3D_CUDA; BP3D_CUDA
cfg = astra_struct('CGLS3D_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = proj_id;


% Create the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);

% Run 150 iterations of the algorithm
% Note that this requires about 750MB of GPU memory, and has a runtime
% in the order of 10 seconds.
astra_mex_algorithm('iterate', alg_id, 150);

% Get the result
rec = astra_mex_data3d('get', rec_id);
figure, imshow(squeeze(rec(:,:,65)),[]);


% Clean up. Note that GPU memory is tied up in the algorithm object,
% and main RAM in the data objects.
astra_mex_algorithm('delete', alg_id);
astra_mex_data3d('delete', rec_id);
astra_mex_data3d('delete', proj_id);

%% root square error
diff = rec - cube;
diff = diff.^2;
err = sum(diff(:))/numel(rec);
disp(['Error ' num2str(err)]);