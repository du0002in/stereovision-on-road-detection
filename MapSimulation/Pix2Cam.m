%Cam intrinsic matrix
%[x_pix;y_pix;1]=Pix2Cam*[x_cam;y_cam;z_cam;1];
function [Tpix2cam_left, Tpix2cam_right]=Pix2Cam(focal,cc_x_left,cc_x_right,cc_y)
Tpix2cam_left=[-focal 0 cc_x_left;0 -focal cc_y;0 0 1];
Tpix2cam_right=[-focal 0 cc_x_right;0 -focal cc_y;0 0 1];