function rec = fbp_gpu(sinogram, vol_geom, proj_geom, filter_type)

% Set default filter type to Ram-Lak if unspecified
if nargin < 4 || isempty(filter_type)
    filter_type = 'Ram-Lak';
end

% Create a data object for the reconstruction
rec_id = astra_mex_data2d('create', '-vol', vol_geom);

% Create a data object to hold the sinogram data
sinogram_id = astra_mex_data2d('create', '-sino', proj_geom, sinogram);

% Create configuration for reconstruction algorithm
cfg = astra_struct('FBP_CUDA');
cfg.ReconstructionDataId = rec_id;
cfg.ProjectionDataId = sinogram_id;

% possible values for FilterType:
% none, ram-lak, shepp-logan, cosine, hamming, hann, tukey, lanczos,
% triangular, gaussian, barlett-hann, blackman, nuttall, blackman-harris,
% blackman-nuttall, flat-top, kaiser, parzen
cfg.FilterType = filter_type;

% Create and run the algorithm object from the configuration structure
alg_id = astra_mex_algorithm('create', cfg);
astra_mex_algorithm('run', alg_id);

% Get the result
rec = astra_mex_data2d('get', rec_id);

% Fix inconsistent scaling and orientation of GPU reconstruction
cor_fac = proj_geom.DistanceOriginSource * ...
    (proj_geom.DistanceOriginSource + proj_geom.DistanceOriginDetector);
rec = cor_fac*rot90(rec,2);

% Clean up. 
astra_mex_algorithm('delete', alg_id);
astra_mex_data2d('delete', rec_id);
astra_mex_data2d('delete', sinogram_id);
