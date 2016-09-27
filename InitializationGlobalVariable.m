%For Kalman Filter
%system performance variance
var_v=0.03;
var_phai=0.001;
%measurement variance
var_w_x=0.0005;
var_w_y=0.0005;
var_w_theta=0.0003;
var_w_v=0.01;
var_w_phai=0.001;

% % % % %For camera intrinsic parameters
% % % % focal=717;
% % % % ox=307;
% % % % oy=257;
% % % % Y=1.547; %Previouse is 1.525
% % % % phi=-47.85/180*pi;

%For car dimentions
%The real minimum turning radius is 2.85m and corresponding to a turning angle of 28.5 degree.
%By assigning 3.1m turnning radius, there will be more buffers in the
%turning angle to turn the vehicle close to the path
Rmin=3.6; 
W_car=0.995;
L_car=1.94;
L_wheels=1.28;
L_wheel_cam=0.675; %horizontal distance between stero camera and the rear wheel

% % % % %For car park slot dimentions
% % % % W_slot=1.4; %inner width
% % % % L_slot=2.55; %outer length
% % % % L_support_points=0;
% % % % W_support_points=0;

% path=[];
P_No=300;
connected_pixel=4.5;