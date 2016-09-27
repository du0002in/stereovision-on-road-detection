function [xc, d_theta, k, z_tangent]=solve_xc_d_theta(abc,x0,z0,theta0)
a=abc(1);b=abc(2);c=abc(3);
co=[2*a*a, 3*a*b, b*b+2*a*(c-x0)+1, b*c-b*x0-z0];
r=roots(co);
r=r(r==real(r));
rl=r;
d_sqr=(a*rl.^2+b*rl+c-x0).^2+(rl-z0).^2;
[~,ind]=min(d_sqr);
% xc=sqrt(d_sqr_min);
z_tangent=rl(ind);
x_tangent=a*z_tangent^2+b*z_tangent+c;
d_theta=atan(1/(2*a*z_tangent+b));
if d_theta<0
    d_theta=d_theta+pi;
end
%The following definition is from the paper on reverse parking
xc=-sin(d_theta)*(x0-x_tangent)+cos(d_theta)*(z0-z_tangent);
d_theta=-d_theta+theta0;
% d_theta=-d_theta;

% k=2*a/(1+(2*a*z_tangent+b)^2)^1.5;
k=2*a;
% phi_d=atan(L_wheels*k);