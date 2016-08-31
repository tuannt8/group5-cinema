function sinogram_cropped = crop_sinogram(sinogram, crop_size)

sinogram_cropped = sinogram(:,crop_size+1:end-crop_size);