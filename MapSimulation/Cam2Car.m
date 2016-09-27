%Take camera pitch and height as input
%Output the transform from cam to car coordinate or cam extrinsic matrix
%[x_cam;y_cam;z_cam;1]=Tcam2car*[x_car;y_car;z_car;1]
function [Tcam2car_left,Tcam2car_right]=Cam2Car(phi,H,dist_left_right,L_wheel_cam)
%rotate cam cor. around its x axis by -phi
R1=[1 0 0;0 cos(-phi) -sin(-phi);0 sin(-phi) cos(-phi)];
%move the new cor. along its z axis by -L_wheel_cam and its y axis by -H
%for left cam, move along its x axis by -0.5*dist_left_right;
%for right cam, move along its x axis by +0.5*dist_left_right;
T2_left=[-0.5*dist_left_right;-H;-L_wheel_cam];
T2_right=[0.5*dist_left_right;-H;-L_wheel_cam];
Tcam2car_left=[R1 R1*T2_left;0 0 0 1];
Tcam2car_right=[R1 R1*T2_right;0 0 0 1];