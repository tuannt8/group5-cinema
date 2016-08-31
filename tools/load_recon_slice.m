function X = load_recon_slice(filename,dims)

% Read in unsigned 16-bit integers
X = load_bin(filename,'*uint16');

% Fix orientation
X = rot90(reshape(X,dims)',2);