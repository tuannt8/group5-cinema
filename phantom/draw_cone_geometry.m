function draw_cone_geometry(geo, obj_size)


    % Draw cone
    detect_size = geo.DetectorRowCount;
    h = detect_size/2;
    dis_source = geo.DistanceOriginSource;
    dis_detect = geo.DistanceOriginDetector;
    
    figure; hold on;
    plot3([0,-dis_source],[0, 0], [0,0]); 
    plot3([0, dis_detect],[0, 0], [0,0]); 
    
    plot3([-dis_source, dis_detect],[0, h], [0,h]);
    plot3([-dis_source, dis_detect],[0, -h], [0,h]);
    plot3([-dis_source, dis_detect],[0, h], [0,-h]);
    plot3([-dis_source, dis_detect],[0, -h], [0,-h]);
    
    plot3([dis_detect,dis_detect],[-h, -h], [-h,h]);
    plot3([dis_detect,dis_detect],[-h, h], [-h,-h]);
    plot3([dis_detect,dis_detect],[h, -h], [h,h]);
    plot3([dis_detect,dis_detect],[h, h], [h,-h]);
    
    x=[0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*obj_size;
    y=[0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*obj_size;
    z=[0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*obj_size;
    
    x = x - obj_size/2;
    y = y - obj_size/2;
    z = z - obj_size/2;
    
for i=1:6
    h=patch(x(:,i),y(:,i),z(:,i),'k');
    set(h,'edgecolor','w')
end
    
    xlim([-dis_source*1.3, dis_detect*2]);
    axis equal;
end