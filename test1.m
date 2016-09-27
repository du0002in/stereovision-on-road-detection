load 'odemotry4.mat';
ite_s=2066; ite_e=2153;
x(ite_s)=0;y(ite_s)=0;theta(ite_s)=pi/2;
v(ite_s)=act_v(ite_s);
phai(ite_s)=act_phai(ite_s);
xdot(ite_s)=cos(theta(ite_s))*v(ite_s);
ydot(ite_s)=sin(theta(ite_s))*v(ite_s);
thetadot(ite_s)=-tan(phai(ite_s))/1.28*v(ite_s);
for ii=(ite_s+1):ite_e
    t=act_t(ii);
    x(ii)=x(ii-1)+xdot(ii-1)*t;
    y(ii)=y(ii-1)+ydot(ii-1)*t;
    theta(ii)=theta(ii-1)+thetadot(ii-1)*t;
    xdot(ii)=cos(theta(ii-1))*v(ii-1);
    ydot(ii)=sin(theta(ii-1))*v(ii-1);
    thetadot(ii)=-tan(phai(ii-1))/1.28*v(ii-1);
    v(ii)=act_v(ii-1);
    phai(ii)=act_phai(ii-1);
end