function sinogram = extract_sinogram_bin(pathname, fileformat, row_num, dims)

%% Get list of files and number of files
filelist = dir(fullfile(pathname,fileformat));
num_files = length(filelist);
num_detector_cols = dims(2);

%% Load first file to get dimensions for initializing sinogram
%X = imread(fullfile(pathname,filelist(1).name));
%num_detectors = size(X,2);

%% Initialize empty sinogram and store selected row of first image.
sinogram = zeros(num_files,num_detector_cols);
%sinogram(1,:) = X(row_num,:);

%% Read all remaining files, extract row and store in sinogram
for k = 1:num_files
    X = load_projection(fullfile(pathname,filelist(k).name),dims);
    %X = imread(fullfile(pathname,filelist(k).name));
    sinogram(k,:) = X(row_num,:);
    if ~mod(k,100)
        fprintf('Loaded file %d of %d.\n',k,num_files); % Display progress
    end
end

%% Values are in percentage transmission I/I_0, convert to [0,1].
sinogram = 0.01*sinogram;

fprintf('DONE: Successfully loaded %d files.\n',num_files); % Display progress
