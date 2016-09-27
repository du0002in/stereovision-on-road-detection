function Particle=p_initialization
Particle=zeros(8,P_No);
Particle(3,:)=pi/2; %theta is 90 degree
[v, phai]=get_measurements;
if v~=0
    v0=v+sqrt(var_v).*randn(1,P_No);
else
    v0=zeros(1,P_No);
end
phai0=phai+sqrt(var_phai).*randn(1,P_No);
Particle(4,:)=cos(Particle(3,:)).*v0;
Particle(5,:)=sin(Particle(3,:)).*v0;
Particle(6,:)=tan(phai0).*v0/L_wheels;
Particle(7,:)=v0;
Particle(8,:)=phai0;