function [Xw_c, Zw_c]=conjugate_points(pre_Xw, pre_Zw, lane_width, leftrightlane, pre_leftrightlane)
if leftrightlane==pre_leftrightlane
    Xw_c=pre_Xw; Zw_c=pre_Zw;
    return;
else
    A=[pre_Zw.^2; pre_Zw; ones(1,length(pre_Zw))]';
    B=pre_Xw';
    abc=pinv(A)*B;
    slope_x=2*abc(1)*pre_Zw+abc(2);
    theta=atan(slope_x);
    if leftrightlane==1
        Xw_c=pre_Xw-lane_width*cos(theta);
        Zw_c=pre_Zw+lane_width*sin(theta);
    else
        Xw_c=pre_Xw+lane_width*cos(theta);
        Zw_c=pre_Zw-lane_width*sin(theta);
    end
end
end