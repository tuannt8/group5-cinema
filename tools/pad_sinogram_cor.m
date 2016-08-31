function sinogram_cor = pad_sinogram_cor(sinogram, cor_padsize)

if cor_padsize > 0
    sinogram_cor = [repmat(sinogram(:,1),1,cor_padsize), sinogram];
else
    sinogram_cor = [sinogram, repmat(sinogram(:,end),1,-cor_padsize)];
end