load 'map_cordinate2';
Lane_marking_width_boundary=0.15;
Lane_marking_width_center=0.15;
Lane_width=3.5;
Lane_marking_length=2.0;
Lane_marking_gap=4.0;
distmat=[Lane_width+0.5*Lane_marking_width_boundary, ...
        Lane_width-0.5*Lane_marking_width_boundary, ...
        0.5*Lane_marking_width_center];
distmat_num=length(distmat);
[Lane_marking_point_x(1,:), Lane_marking_point_y(1,:)]=getlanemarkingpoint(map_x([end,1,2]), ...
        map_y([end,1,2]),distmat);
for i=2:(length(map_x)-1)
    [xtmp, ytmp]=getlanemarkingpoint(map_x(i-1:i+1),map_y(i-1:i+1), ...
        [Lane_width+0.5*Lane_marking_width_boundary, ...
        Lane_width-0.5*Lane_marking_width_boundary, ...
        0.5*Lane_marking_width_center]);
    d1=(xtmp(1)-Lane_marking_point_x(i-1,1))^2+(ytmp(1)-Lane_marking_point_y(i-1,1))^2;
    d2=(xtmp(1)-Lane_marking_point_x(i-1,distmat_num+1))^2+(ytmp(1)-Lane_marking_point_y(i-1,distmat_num+1))^2;
    if d1<d2
        Lane_marking_point_x(i,:)=xtmp;
        Lane_marking_point_y(i,:)=ytmp;
    else
        Lane_marking_point_x(i,1:distmat_num)=xtmp(distmat_num+1:end);
        Lane_marking_point_x(i,distmat_num+1:2*distmat_num)=xtmp(1:distmat_num);
        Lane_marking_point_y(i,1:distmat_num)=ytmp(distmat_num+1:end);
        Lane_marking_point_y(i,distmat_num+1:2*distmat_num)=ytmp(1:distmat_num);
    end
end
[Lane_marking_point_x(length(map_x),:), Lane_marking_point_y(length(map_x),:)]=getlanemarkingpoint(map_x([1,end,end-1]),map_y([1,end,end-1]), ...
        [Lane_width+0.5*Lane_marking_width_boundary, ...
        Lane_width-0.5*Lane_marking_width_boundary, ...
        0.5*Lane_marking_width_center]);

    
    
    %the lane marking map generated is not perfect from this code. some
    %couter occurs at sharp turning points.