function inlier=simpleransac(x,y)
inliernum=0;
for i=1:1000
    ind=randi(length(x),1,2);
    u=x(ind);
    v=y(ind);
    c=-(v(1)-v(2))/(u(1)-u(2));
    fp=-c*u(1)-v(1);
    deno=sqrt(c^2+1);
    ds=abs(c*x+y+fp)/deno;
    inliertmp=(ds<2.0);
    if sum(inliertmp)>inliernum
        inliernum=sum(inliertmp);
        inlier=inliertmp;
    end
end