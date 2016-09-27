% inverse perspective mapping
% img is the retified and half resoluted image, p is view angle Phi, H is cam height, focal is
% cam focal length, xxc and yyc are cam coordinate grid defined as
% following:
% [ximg,yimg]=meshgrid(1:640,1:2:480);
% yyc=-yimg+cc_y;
% xxc=-ximg+cc_x_left;
function ipm_img=ipm(img,p,H,focal,xxc,yyc)
tic;
Zw=(yyc*H*sin(p)+focal*H*cos(p))./(focal*sin(p)-yyc*cos(p));
Xw=xxc.*(H*sin(p)+Zw*cos(p))./focal;
toc
ind=((abs(Xw)<3) & (Zw>1) & (Zw<10));
ind=find(ind==1);
ipm_img=zeros(225,300);
row_ind=round(-(Zw(ind)-1)*25+225);
row_ind(row_ind<1)=1; row_ind(row_ind>225)=225;
col_ind=round(-(Xw(ind)--3)*50+300);
col_ind(col_ind<1)=1; col_ind(col_ind>300)=300;
ipm_ind=sub2ind([225 300],row_ind,col_ind);
ipm_img(ipm_ind)=img(ind);
toc
end