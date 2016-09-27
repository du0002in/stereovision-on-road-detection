function [xtmpl2,ytmpl2,xtmpr2,ytmpr2,num_in,dumy2]=fit_other_lane_zebra(leftrightlane,dumy1,xtmpl,ytmpl,xtmpr,ytmpr,D_pre,H_pre,ximg,yimg)
global_varibles;
ImgSize2=640;ImgSize1=480;
xtmpl2=[];ytmpl2=[];xtmpr2=[];ytmpr2=[];num_in=0;
if leftrightlane==0
    li_lb=ones(1,length(ytmpl));
    li_rb=xtmpl+round(46.8-2.26*ytmpl);
    ri_lb=ones(1,length(ytmpr));
    ri_rb=xtmpr+round(46.8-2.26*ytmpr);
else
    li_lb=xtmpl-round(46.8-2.26*ytmpl);
    li_rb=ImgSize2*ones(1,length(ytmpl));
    ri_lb=xtmpr-round(46.8-2.26*ytmpr);
    ri_rb=ImgSize2*ones(1,length(ytmpr));
end
dumy2=C_in_region(dumy1,li_lb,li_rb,ri_lb+ImgSize2+30,ri_rb+ImgSize2+30,-ytmpr(1),240);
dumy2=dumy1-dumy2;

dumy2_l=dumy2(:,1:ImgSize2);
dumy2_r=dumy2(:,(end-ImgSize2+1):end);
[yl, xl]=find(dumy2_l==1);
yl=-(2*yl-1)+cc_y;
xl=-xl+cc_x_left;
[yr, xr]=find(dumy2_r==1);
yr=-(2*yr-1)+cc_y;
xr=-xr+cc_x_right;

Datapointl=[yl xl]';
Datapointr=[yr xr]';
% if isempty(Datapointl) || isempty(Datapointr)
if length(Datapointl)<4 || length(Datapointr)<4
    num_in=0;
else
    [linel, liner, para, y_int]=P_Zebra_RANSAC_new(Datapointl, Datapointl, Datapointr, Datapointr, D_pre, H_pre, 3000, focal, dist_left_right, cc_y);
    num_in=sum(linel)+sum(liner);
end

if (num_in>50) && (sum(yl(linel==1)>y_int(1))+sum(yr(liner==1)>y_int(2)))>30 && ...
            (sum(yl(linel==1)<=y_int(1))+sum(yr(liner==1)<=y_int(2)))>30
    y_int_ave=mean(y_int);
    top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
    btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even
    ytmpl2=top_row:-2:btm_row;
    xtmpl2=(-para(3,2)-2*para(2,2)*ytmpl2(ytmpl2>=y_int_ave))/(2*para(1,2));
    xtmpl2=[xtmpl2 (-para(3,1)-2*para(2,1)*ytmpl2(ytmpl2<y_int_ave))/(2*para(1,1))];
    ytmpr2=top_row:-2:btm_row;
    xtmpr2=(-para(5,2)-2*para(4,2)*ytmpr2(ytmpr2>=y_int_ave))/(2*para(1,2));
    xtmpr2=[xtmpr2 (-para(5,1)-2*para(4,1)*ytmpr2(ytmpr2<y_int_ave))/(2*para(1,1))];
    ytmpl2=-(-ytmpl2+cc_y+1)/2;
    xtmpl2=-(xtmpl2-cc_x_left);
    ytmpr2=-(-ytmpr2+cc_y+1)/2;
    xtmpr2=-(xtmpr2-cc_x_right);
    
    ytmpl2(round(xtmpl2)<=0)=[];
    xtmpl2(round(xtmpl2)<=0)=[];
    ytmpl2(round(xtmpl2)>ImgSize2)=[];
    xtmpl2(round(xtmpl2)>ImgSize2)=[];
    
    ytmpr2(round(xtmpr2)<=0)=[];
    xtmpr2(round(xtmpr2)<=0)=[];
    ytmpr2(round(xtmpr2)>ImgSize2)=[];
    xtmpr2(round(xtmpr2)>ImgSize2)=[];
    
    if ~isempty(ximg)
        ytest=-ytmpr2; xtest=xtmpr2;
        ximg1=ximg(1:4); yimg1=yimg(1:4)/2;
        min_x=max([1, min(ximg1)]);   max_x=min([ImgSize2, max(ximg1)]);
        min_y=max([1, min(yimg1)]);   max_y=min([ImgSize1/2, max(yimg1)]);
        indtmp=(xtest>=min_x & xtest<=max_x & ytest>=min_y & ytest<=max_y);
        if sum(indtmp)~=0
            xtest=xtest(indtmp==1); ytest=ytest(indtmp==1);
            if C_pass_letter(xtest,ytest,ximg1,yimg1)==1
                xtmpl2=[];ytmpl2=[];xtmpr2=[];ytmpr2=[];
                return;
            end
        end
    end 
else
    num_in=0;
end