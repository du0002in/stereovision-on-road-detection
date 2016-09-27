function [marking_x, marking_y]=getlanemarkingpoint(x,y,dist)
if x(1)==x(3)
    e=0; c=1; f=-2*x(1); fp=f/2;
else
    c=-(y(1)-y(3))/(x(1)-x(3));
    f=2*(-c*x(1)-y(1));
    e=1;fp=f/2;
end
v_c=-e; v_e=c; v_fp=-v_c*x(2)-v_e*y(2);
if v_e~=0
    ax=1+(v_c/v_e)^2;
    bx=2*((v_fp/v_e+y(2))*v_c/v_e-x(2));
    cx=(v_fp/v_e+y(2))^2+x(2)^2-dist.^2;
    xtmp=(-bx+sqrt(bx^2-4*ax*cx(1)))/(2*ax);
    ytmp=-v_c/v_e*xtmp-v_fp/v_e;
    dtmp=c*xtmp+e*ytmp+fp;
%     if ((-v_c/v_e)<0 && dtmp>0) || ((-v_c/v_e)>=0 && dtmp<=0)
        marking_x(1:length(dist))=(-bx+sqrt(bx^2-4*ax*cx))/(2*ax);
        marking_x((length(dist)+1):(2*length(dist)))=(-bx-sqrt(bx^2-4*ax*cx))/(2*ax);
%     else
%         marking_x(1:length(dist))=(-bx-sqrt(bx^2-4*ax*cx))/(2*ax);
%         marking_x((length(dist)+1):(2*length(dist)))=(-bx+sqrt(bx^2-4*ax*cx))/(2*ax);
%     end
    marking_y=-v_c/v_e*marking_x-v_fp/v_e;
    
%     
%     if (-v_c/v_e)>0
%         marking_x(1:length(dist))=(-bx+sqrt(bx^2-4*ax*cx))/(2*ax);
%         marking_x((length(dist)+1):(2*length(dist)))=(-bx-sqrt(bx^2-4*ax*cx))/(2*ax);
%     else
%         marking_x(1:length(dist))=(-bx-sqrt(bx^2-4*ax*cx))/(2*ax);
%         marking_x((length(dist)+1):(2*length(dist)))=(-bx+sqrt(bx^2-4*ax*cx))/(2*ax);
%     end
%     marking_y=-v_c/v_e*marking_x-v_fp/v_e;
else
    ay=1;
    by=-2*y(2);
    cy=(v_fp/v_c+x(2))^2+y(2)^2-dist.^2;
    marking_x(1:(2*length(dist)))=-v_fp/v_c*ones([1,(2*length(dist))]);
    marking_y(1:length(dist))=(-by+sqrt(by^2-4*ay*cy))/(2*ay);
    marking_y((length(dist)+1):(2*length(dist)))=(-by-sqrt(by^2-4*ay*cy))/(2*ay);
end