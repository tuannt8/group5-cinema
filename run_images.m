clear;
close all;

load('fdk_3d_11slices_center_508_subsampled3.mat');


center_slice = 508;
slice_range = (-5:1:5) + center_slice;

figure
for i = 1:size(reconstruction,3)
    show_image(rot90(squeeze(reconstruction(:,:,i)),-1));
    title(['FDK, Slice = ' num2str(slice_range(i))]);
    drawnow
    pause(1)
end

load('cgls3d_10iter_11slices_center_508.mat');

center_slice = 512;
slice_range = (-5:1:5) + center_slice;


figure
for i = 1:size(reconstruction_cgls,3)
    show_image(rot90(squeeze(reconstruction_cgls(:,:,i)),-1))
    title(['CGLS, 10 iters, Slice = ' num2str(slice_range(i))]);
    drawnow
    pause(1)
end

load('cgls3d_20iter_11slices_center_508.mat');

center_slice = 512;
slice_range = (-5:1:5) + center_slice;


figure
for i = 1:size(reconstruction_cgls2,3)
    show_image(rot90(squeeze(reconstruction_cgls2(:,:,i)),-1))
    title(['CGLS, 20 iters, Slice = ' num2str(slice_range(i))]);
    drawnow
    pause(1)
end