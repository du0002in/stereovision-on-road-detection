function [xc, d_theta]=bench_mark_algo1(img,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y)

[ipm,~]=C_ipm_no_rotate2(img,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);

[r,c]=find(ipm==1);
Xw=(150-c)/50;
Zw=(500-r)/50;
if length(Xw)>=5
[~,para]=C_fit_single_line_Threshold(Zw,Xw,0.04);
abc=-[0; para(2); para(3)]/para(1);
[xc, d_theta]=solve_xc_d_theta(abc,0,0,pi/2);
else
    xc=nan; d_theta=nan;
end
end