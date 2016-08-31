function [] = plot_roi(N,num_detector_pixels)

theta = linspace(0,2*pi,1000);
x = 0.5 + 0.5*N + 0.5*num_detector_pixels*cos(theta);
y = 0.5 + 0.5*N + 0.5*num_detector_pixels*sin(theta);

plot(x,y,'-r')