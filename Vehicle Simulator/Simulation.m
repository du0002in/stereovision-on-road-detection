clear all;
sendtarget(0,0);
global_varibles;
InitializationGlobalVariable;
P_s=0.2; Q_s=0.5;

%car coordinate:    ini_xc, ini_d_theta, plotCar, abc, VehicleSimulator,
%inner_thread, act_*, mea_d_theta, mea_xc,
%global coordin:    z_dumy, x_dumy, plot, 

ini_xc=-0.8057; ini_d_theta=-0.0863;
z_dumy=[0, 60]; x_dumy=-(tan(ini_d_theta)*z_dumy+ini_xc/cos(ini_d_theta));
figure;plot(x_dumy,z_dumy); hold on; plotCar(pi/2,0,0,0); hold off;axis equal;
abc=[0 tan(ini_d_theta) ini_xc/cos(ini_d_theta)];
VehicleSimulator('Ini',ini_xc, ini_d_theta);pause(1);
inner_thread('Ini',ini_xc,ini_d_theta,abc,L_wheels,Q_s,P_s,0,2.0);
for i=1:50
    pause(0.45);
    [act_x(i), act_z(i), act_theta(i), act_phi(i), act_xc(i), act_d_theta(i)]=getActPose;
    plot(x_dumy,z_dumy); hold on;plot(-act_x,act_z,'r'); plotCar(pi-act_theta(i),-act_x(i),act_z(i),act_phi(i)); hold off;axis equal;
%     mea_d_theta(i)=act_d_theta(i)+normrnd(0,0.0334,[1 1]);
%     mea_xc(i)=act_xc(i)+normrnd(0,0.0914,[1 1]);
    mea_d_theta(i)=act_d_theta(i);
    mea_xc(i)=act_xc(i);
    abc=[0 tan(mea_d_theta(i)) mea_xc(i)/cos(mea_d_theta(i))];
    inner_thread('Con',mea_xc(i), mea_d_theta(i),abc,L_wheels,Q_s,P_s,0,2.0);    
end

VehicleSimulator('End', -5, 0);
inner_thread('End',-3.0,0,[0 0 -3.0],L_wheels,0.5,0.2,0,0.6);
hold on; plotCar(pi/2,0,0,0); hold off;
hold on;plot(x_dumy-3.5, z_dumy,'b'); hold off;