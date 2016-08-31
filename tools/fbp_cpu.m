function rec = fbp_cpu(sinogram, proj_id)

% Get projection geometry from projector
proj_geom = astra_mex_projector('projection_geometry', proj_id);
vol_geom = astra_mex_projector('volume_geometry', proj_id);

% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);

% Create a data object to hold the sinogram data
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);

% Create configuration for reconstruction algorithm
cfg = astra_struct('FBP');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;
cfg.ProjectorId = proj_id;

% Create and run the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', alg_id);

% Get the result
rec = astra_mex_data2d('get', rec_id);

% Clean up. 
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
