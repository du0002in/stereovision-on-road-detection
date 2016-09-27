function [Xw_new, Zw_new]=best_line(Xw,Zw,Xc,D_theta,leftrightlane,pre_leftrightlane,lane_width, num_in, Xw1, Zw1)
var_w_xc=(450/min([num_in,450])*0.0914)^2; 
var_w_d_theta=(450/min([num_in,450])*0.0334).^2;
if leftrightlane~=pre_leftrightlane
    if pre_leftrightlane==0
        Xc=Xc-lane_width;
    elseif pre_leftrightlane==1
        Xc=Xc+lane_width;
    end
end

A=[Zw; ones(1, length(Zw))]';
B=Xw';
bc=pinv(A)*B;
abc=[0; bc];
[xc_m, d_theta_m,~]=solve_xc_d_theta(abc,0,0,pi/2);
w(1)=exp(-0.5*(1.0*(xc_m-Xc).^2/var_w_xc+1.0*(d_theta_m-(-0.0133)-D_theta).^2/var_w_d_theta));

% ctr_x=0.5*(Xw(1)+Xw(end)); ctr_z=0.5*(Zw(1)+Zw(end));
b=bc(1); c=bc(2);
ctr_z=-(b*(c-Xw1)-Zw1+sqrt((b*(c-Xw1)-Zw1)^2-(b*b+1)*((c-Xw1)^2+Zw1^2-4.0)))/(b*b+1);
ctr_x=b*ctr_z+c;

b1=bc(1); b1=tan(atan(b1)+0.1);
c1=ctr_x-b1*ctr_z;
[xc_m, d_theta_m,~]=solve_xc_d_theta([0;b1;c1],0,0,pi/2);
w(2)=exp(-0.5*(1.0*(xc_m-Xc).^2/var_w_xc+1.0*(d_theta_m-(-0.0133)-D_theta).^2/var_w_d_theta));

b=bc(1); b=tan(atan(b)-0.1);
c=ctr_x-b*ctr_z;
[xc_m, d_theta_m,~]=solve_xc_d_theta([0;b;c],0,0,pi/2);
w(3)=exp(-0.5*(1.0*(xc_m-Xc).^2/var_w_xc+1.0*(d_theta_m-(-0.0133)-D_theta).^2/var_w_d_theta));

[~,ind]=max(w);
if ind==1
    Xw_new=Xw; Zw_new=Zw;
elseif ind==2
    Zw_new=Zw;
    Xw_new=b1*Zw_new+c1;
else
    Zw_new=Zw;
    Xw_new=b*Zw_new+c;
end

end