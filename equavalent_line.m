function [X, Y, Z, abc]=equavalent_line(break_location, Xw, Zw)
Atmp=[Zw(1:(break_location-1)); ones(1,length(Zw(1:(break_location-1))))]';
Btmp=Xw(1:(break_location-1))';
bc_top=pinv(Atmp)*Btmp;
% angle_top=atan(bc_top(1));

Atmp=[Zw(break_location:end); ones(1,length(Zw(break_location:end)))]';
Btmp=Xw(break_location:end)';
bc_btm=pinv(Atmp)*Btmp;
% angle_btm=atan(bc_btm(1));

zw_int=(bc_top(2)-bc_btm(2))/(bc_btm(1)-bc_top(1));
xw_int=bc_top(1)*zw_int+bc_top(2);

b=bc_top(1); c=bc_top(2);
ztop=(-(b*(c-xw_int)-zw_int)+sqrt((b*(c-xw_int)-zw_int)^2-(b^2+1)*(zw_int^2+(c-xw_int)^2-2.05^2)))/(b^2+1);
if ztop<zw_int
    ztop=(-(b*(c-xw_int)-zw_int)-sqrt((b*(c-xw_int)-zw_int)^2-(b^2+1)*(zw_int^2+(c-xw_int)^2-2.05^2)))/(b^2+1);
end
xtop=b*ztop+c;

b=bc_btm(1); c=bc_btm(2);
zbtm=(-(b*(c-xw_int)-zw_int)-sqrt((b*(c-xw_int)-zw_int)^2-(b^2+1)*(zw_int^2+(c-xw_int)^2-2.05^2)))/(b^2+1);
if zbtm>zw_int
    zbtm=(-(b*(c-xw_int)-zw_int)+sqrt((b*(c-xw_int)-zw_int)^2-(b^2+1)*(zw_int^2+(c-xw_int)^2-2.05^2)))/(b^2+1);
end
xbtm=b*zbtm+c;

Atmp=[ztop zbtm;1 1]';
Btmp=[xtop xbtm]';
bc=pinv(Atmp)*Btmp;
abc=[0;bc];

Z=Zw; Y=zeros(1,length(Zw));
X=abc(2)*Z+abc(3);
end