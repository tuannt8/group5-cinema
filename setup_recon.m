
%% Add tools subdirectory to MATLAB path
addpath tools

%% Of the two gpus available, choose a random one based on system time
cc = clock;
gpu_index = mod(round(cc(6)),2)+1
gpuDevice(gpu_index)

%% Add ASTRA to MATLAB path if not already present
