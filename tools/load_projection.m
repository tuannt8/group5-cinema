function X = load_projection(filename,dims)

X = load_bin(filename,'*float');
X = reshape(X,dims)';