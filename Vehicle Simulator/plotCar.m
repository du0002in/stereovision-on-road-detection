function plotCar(theta,x,z, phai)
global_varibles;
theta=-pi/2-theta;
R=[cos(theta) 0 sin(theta);...
    0 1 0;...
    -sin(theta) 0 cos(theta)];

%plot car outline
x_car=[0.5*W_car -0.5*W_car -0.5*W_car 0.5*W_car];
z_car=[(L_car-L_wheels)*0.5 (L_car-L_wheels)*0.5 -L_car+(L_car-L_wheels)*0.5 -L_car+(L_car-L_wheels)*0.5];
car_global=R*[x_car; zeros(1,4); z_car];
x_car_global=car_global(1,:)+x;
z_car_global=car_global(3,:)+z;hold on;
% plot([x_car_global x_car_global(1)], [z_car_global z_car_global(1)],'LineWidth',2);
patch(x_car_global, z_car_global,'c');

%plot car front direction
x_d=[-0.5*W_car 0 0.5*W_car];
z_d=[0, -L_car+(L_car-L_wheels)*0.5, 0];
d_global=R*[x_d; zeros(1,3); z_d];
x_d_global=d_global(1,:)+x;
z_d_global=d_global(3,:)+z;
plot([x_d_global x_d_global(1)],[z_d_global z_d_global(1)]);
x_d=[0 0];z_d=[0, -L_car+(L_car-L_wheels)*0.5];
d_global=R*[x_d; zeros(1,2); z_d];
x_d_global=d_global(1,:)+x;
z_d_global=d_global(3,:)+z;
plot([x_d_global],[z_d_global]);
%plot car wheel steering
Rear_right_x=[W_car*0.5-0.06 W_car*0.5-0.06];
Rear_right_z=[0.25 -0.25]; %0.25 is the radius of the wheel
Rear_left_x=-[W_car*0.5-0.06 W_car*0.5-0.06];
Rear_left_z=[0.25 -0.25];
Front_right_x=[W_car*0.5-0.06-.25*sin(phai) W_car*0.5-0.06+0.25*sin(phai)];
Front_right_z=[-L_wheels+0.25*cos(phai) -L_wheels-0.25*cos(phai)];
Front_left_x=[-W_car*0.5+0.06-.25*sin(phai) -W_car*0.5+0.06+0.25*sin(phai)];
Front_left_z=[-L_wheels+0.25*cos(phai) -L_wheels-0.25*cos(phai)];
wheel_matrix_x=[Rear_right_x Rear_left_x Front_right_x Front_left_x];
wheel_matrix_z=[Rear_right_z Rear_left_z Front_right_z Front_left_z];
wheel_global=R*[wheel_matrix_x; zeros(1,8); wheel_matrix_z];
wheel_global_x=wheel_global(1,:)+x;
wheel_global_z=wheel_global(3,:)+z;
plot(wheel_global_x([1,2]),wheel_global_z([1,2]),'k','LineWidth',4);
plot(wheel_global_x([3,4]),wheel_global_z([3,4]),'k','LineWidth',4);
plot(wheel_global_x([5,6]),wheel_global_z([5,6]),'k','LineWidth',4);
plot(wheel_global_x([7,8]),wheel_global_z([7,8]),'k','LineWidth',4);

% % % % %plot vision range
% % % % x_vision=[0.68 1.10 -1.10 -0.68];
% % % % z_vision=[0.7 2.8 2.8 0.7]+L_wheel_cam;
% % % % vision_global=R*[x_vision; zeros(1,4); z_vision];
% % % % x_vision_global=vision_global(1,:)+x;
% % % % z_vision_global=vision_global(3,:)+z;
% % % % plot(x_vision_global, z_vision_global,'b--');
% % % % hold off;