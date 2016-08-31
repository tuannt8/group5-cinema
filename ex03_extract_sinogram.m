%% Beware: Clear existing variables
clear

%% Parameters to specify

% Case FoV2_5mm or LFoV
case_id = 'FoV2_5mm';

% Which row of projection images to extract. The central slice is 508 for
% FoV2_5mm and 512 for LFoV
row_num = 1000;

% Projections to show
proj_idx = [1,100,400,800,1200, 2400, 4000, 4600];

% Whether to save generated sinogram and savefilename
do_save = true;

%% Paths to data, file format, name of file for saving sinogram

% Path to data
pathname_data = sprintf(...
    '/work1/jakj/cinemax2/data/CaseGFRP/%s/Projections/bin/',case_id);

% Geometry file
geometry_file = sprintf(...
    '/work1/jakj/cinemax2/data/CaseGFRP/%s/geometry_%s.txt',...
    case_id,case_id);

% File format
file_format = '*.bin';

% Name of file to save sinogram to
savefile = fullfile('coursesinos',...
    sprintf('%s_sino%04d.mat',case_id,row_num));

%% Extract only projection dimensions from geometry file, needed to load
[~, ~, ~, dims] = load_geometry(geometry_file);

%% Load and display a few projections as images

% Get list of bin-files
filelist = dir(fullfile(pathname_data,file_format));

% Load specified projections and display
Xall = zeros(dims(1),dims(2),length(proj_idx));
for k = 1:length(proj_idx)
    Xall(:,:,k) = load_projection(fullfile(pathname_data,...
        filelist(proj_idx(k)).name), dims);
    figure
    show_image(Xall(:,:,k))
    title(sprintf('Projection %d', proj_idx(k)))
    colorbar
    hold on
    plot(xlim,row_num*[1,1],'-r')
end

%% Plot the rows extracted for sinogram as profiles
figure
plot(squeeze(Xall(row_num,:,:)))
legend(num2str(proj_idx'),'location','southeast')
drawnow

%% Load same single row of all tif files in the specified path

% Read from binary files with projection images
sinogram = extract_sinogram_bin(pathname_data, file_format, row_num, dims);

%% Display extracted sinogram

figure
imagesc(sinogram)
axis image
colorbar
colormap gray
ylabel('Projection number')
xlabel('Detector pixel')

%% If asked for, save to file.
if do_save
    save(savefile,'sinogram');
end

