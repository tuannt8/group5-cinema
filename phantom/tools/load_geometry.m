function [source_origin, detector_origin, pixel_size, proj_dims, ...
    proj_num, recon_slice_dims, recon_slice_num] = load_geometry(filename)

geometry = importdata(filename);

source_origin = geometry.data(1);
detector_origin = geometry.data(2);
pixel_size = geometry.data(3);

proj_dims = [geometry.data(4),geometry.data(5)];
proj_num = geometry.data(6);

recon_slice_dims = [geometry.data(7),geometry.data(8)];
recon_slice_num = geometry.data(9);