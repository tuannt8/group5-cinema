function X = load_bin(filename,precision)

fid = fopen(filename,'r');
X = fread(fid,inf,precision);
fclose(fid);