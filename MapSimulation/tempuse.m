load 'map_cordinate3';
ori_map_x=map_x; ori_map_y=map_y;
register_c=length(map_x);ctr=0;
stopflag=1;
j=0;
while stopflag==1
    ctr=ctr+1
j=0;
for i=1:length(map_x)-1
    j=j+1;
    d(i)=sqrt((map_x(i+1)-map_x(i))^2+(map_y(i+1)-map_y(i))^2);
    new_map_x(j)=map_x(i);
    new_map_y(j)=map_y(i);
    new_Lane_marking_point_x(j,:)=Lane_marking_point_x(i,:);
    new_Lane_marking_point_y(j,:)=Lane_marking_point_y(i,:);
    if d(i)>=0.1
        j=j+1;
        xtmp=0.5*(map_x(i+1)+map_x(i));
        ytmp=0.5*(map_y(i+1)+map_y(i));
        new_map_x(j)=xtmp;
        new_map_y(j)=ytmp;
        
        xtmp=0.5*(Lane_marking_point_x(i+1,:)+Lane_marking_point_x(i,:));
        ytmp=0.5*(Lane_marking_point_y(i+1,:)+Lane_marking_point_y(i,:));
        new_Lane_marking_point_x(j,:)=xtmp;
        new_Lane_marking_point_y(j,:)=ytmp;
    end
end
new_map_x(j+1)=map_x(end);
new_map_y(j+1)=map_y(end);
new_Lane_marking_point_x(j+1,:)=Lane_marking_point_x(end,:);
new_Lane_marking_point_y(j+1,:)=Lane_marking_point_y(end,:);

map_x=new_map_x; map_y=new_map_y;
Lane_marking_point_x=new_Lane_marking_point_x;
Lane_marking_point_y=new_Lane_marking_point_y;
% if register_c==j
%     stopflag=0;
% else
%     register_c=length(map_x);
% end
log_j(ctr)=j;
if max(d)<=0.1
    stopflag=0;
end
end