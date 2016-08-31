function sinogram = pad_sinogram_cor_3d(sinogram, cor_padsize)

if cor_padsize > 0
    sinogram = cat(2, repmat(sinogram(:,1,:), 1, cor_padsize, 1), sinogram);
else
    sinogram = cat(2, sinogram, repmat(sinogram(:,end,:), 1, -cor_padsize, 1));
end
end