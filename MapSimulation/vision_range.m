% coordinate follows stereo_img(2)
% function [Xg, Yg, Zg]=vision_range(X,Y)
X=-act_x(ii); Y=act_z(ii); simulated_lane_width=3.5;
ind_x=((map_x>=(X-2*simulated_lane_width))&(map_x<=(X+simulated_lane_width*2)));
ind_y=((map_y>=(Y-2*simulated_lane_width))&(map_y<=(Y+simulated_lane_width*2)));
ind_xy=(ind_x)&(ind_y);
ind_xy=find(ind_xy==1);
map_xtmp=map_x(ind_xy);
map_ytmp=map_y(ind_xy);
[IDX,xc]=knnsearch([map_xtmp' map_ytmp'],[X,Y]);
IDX=ind_xy(IDX);
Lane_marking_point_z=zeros(size(Lane_marking_point_x));
FOV_deapth=60; %the deapth of field of view. 600 data points~=30m
if (IDX+FOV_deapth)<=length(map_x)
    FOV_ind=IDX:(IDX+FOV_deapth);
else
    FOV_ind=[IDX:length(map_x) 1:(FOV_deapth-(length(map_x)-IDX))];
end
Xg=Lane_marking_point_x(FOV_ind,:);
Yg=Lane_marking_point_y(FOV_ind,:);
Zg=Lane_marking_point_z(FOV_ind,:);
