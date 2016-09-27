function [line, bedf]=RANSACtest(Datapoint, SubDatapoint)
tcon_pre=0; tcon=0;
p=[Datapoint(2,:); Datapoint(1,:); ones(1,length(Datapoint))];%3xn
for i=1:1000
    tcon=0;
    Randpoint=randperm(length(SubDatapoint),4);
    v=(SubDatapoint(1,Randpoint))';
    u=(SubDatapoint(2,Randpoint))';
    Vmatrix=[v.*v v u ones(4,1)];
    Umatrix=u.*v;
    bedftmp=Vmatrix\Umatrix;
    B=bedftmp(1);E=bedftmp(2);D=bedftmp(3);F=bedftmp(4);
    C=E+B*D; A=F+C*D;
    if u(1)==u(4)
        e=0;c=1;f=-2*u(1);fp=f/2;
    else
        c=-(v(1)-v(4))/(u(1)-u(4));
        e=1;
        f=2*(-c*u(1)-v(1));
        fp=f/2;
    end
    deno=c*c+e*e;
    Cr=[0 -0.5 D/2;-0.5 B E/2; D/2 E/2 F];
    
    Cp=Cr*p;%3xn
    dstmp=(diag(transpose(p)*Cr*p))';
    ds=(dstmp.*dstmp)/4.0./(Cp(1,:).*Cp(1,:)+Cp(2,:).*Cp(2,:));
    [~,index]=find(ds<5);
    tcon=length(index);
    if tcon>tcon_pre
        tcon_pre=tcon;
        bedf=bedftmp;
        line=index;
    end
end