function inliner_index=Hump_Detection(img,Phai,H,num_ite,focal,theta,cc_y,cc_x_left)
% ImgSize2=size(img,2);
[indvl, indul]=find(img==1);
indvl=-(2*indvl-1)+cc_y;
indul=-indul+cc_x_left;
Datapoint=[indvl,indul]';
SubDatapoint=Datapoint;
[inliner, inliner_index,final_theta,final_xc]= ...
P_Hump_RANSAC(Datapoint,SubDatapoint,Phai,H,num_ite,focal,theta,cc_y,cc_x_left);
end