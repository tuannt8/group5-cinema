function [] = disp_geometry(proj_geom, vol_geom)

% Extract parameters
N = vol_geom.GridRowCount;
source_pos = proj_geom.DistanceOriginSource;
det_pos = proj_geom.DistanceOriginDetector;
det_lim = proj_geom.DetectorCount*proj_geom.DetectorWidth/2;

% Plot source
plot(-source_pos,0,'or','markersize',10)
hold on

% Center of rotation position
plot(0,0,'ok','markersize',10)

% Detector
plot([1,1,1]*det_pos,[det_lim,-det_lim,0],'ob','markersize',10)
plot(det_pos*[1,1],det_lim*[-1,1],'-b')

% Square object box
xo = [1,-1,-1, 1,1]*N/2;
yo = [1, 1,-1,-1,1]*N/2;
plot(xo,yo,'-k')

% Source to detector limits
plot([-source_pos,det_pos],[0,det_lim],'-r')
plot([-source_pos,det_pos],[0,-det_lim],'-r')
plot([-source_pos,det_pos],[0,0],'-r')

% Customize
axis equal
xlim([-source_pos-0.5*det_pos,det_pos+0.5*det_pos])
xlim([-source_pos-N,det_pos+N])
legend('X-ray source & beam',...
    'Object box & center of rotation',...
    'Detector','location','northwest')