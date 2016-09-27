function [sel_xtmpl2,sel_ytmpl2,xtmpl2, ytmpl2, sel_xtmpr2,sel_ytmpr2,xtmpr2, ytmpr2, num_inlier,para2,lane_width,dumy2, xtmpl_imagine, xtmpr_imagine,ytmp_imagine]=fit_other_lane(leftrightlane,dumy,para,ytopr,ybtmr,H,Phai,D_pre, ximg,yimg)
global_varibles;
ImgSize2=640;ImgSize1=480;
xtmpl2=[];ytmpl2=[];xtmpr2=[];ytmpr2=[];para2=[];
sel_xtmpl2=[];sel_ytmpl2=[];sel_xtmpr2=[];sel_ytmpr2=[];lane_width=[];
num_inlier=0;
ybtmr=240;

L=[2.2 4.7]; %Width of the road is assumed to vary from 2.7m to 4.7m
Bp=cos(Phai)/H*L;
if leftrightlane==0 %left lane
    if sum(abs(para(:,2)))==0 %model is hypobola
        newBli=para(2,1)+Bp;
        newBri=para(4,1)+Bp;
        newCli=para(3,1)-D_pre*Bp;
        newCri=para(5,1)-D_pre*Bp;
    else %model is line
        neweli=para(2,2)-Bp;
        neweri=para(4,2)-Bp;
        newfli=para(3,2)+2*D_pre*Bp;
        newfri=para(5,2)+2*D_pre*Bp;
    end
else %right lane
    if sum(abs(para(:,2)))==0 %model is hypobola
        Bli=para(2,1)-Bp;
        Bri=para(4,1)-Bp;
        Cli=para(3,1)+D_pre*Bp;
        Cri=para(5,1)+D_pre*Bp;
    else
        eli=para(2,2)+Bp;
        eri=para(4,2)+Bp;
        fli=para(3,2)-2*D_pre*Bp;
        fri=para(5,2)-2*D_pre*Bp;
    end
end

y=(ytopr:1:ybtmr)*(-2)+cc_y+1; %convert to [ccy:ccy-480]x640
for i=1:length(Bp)
    if leftrightlane==0 %left lane
        if sum(abs(para(:,2)))==0 %model is hypobola
            xtmpli_rlane=(para(1,1)./(y-D_pre)+newBli(i)*y+newCli(i));
            xtmpri_rlane=(para(1,1)./(y-D_pre)+newBri(i)*y+newCri(i));
        else %model is line
            xtmpli_rlane=(-neweli(i)*y-newfli(i)/2);
            xtmpri_rlane=(-neweri(i)*y-newfri(i)/2);
        end
        xtmpli(i,:)=xtmpli_rlane; xtmpri(i,:)=xtmpri_rlane;
    else %right lane
        if sum(abs(para(:,2)))==0 %model is hypobola
            xtmpli_llane=(para(1,1)./(y-D_pre)+Bli(i)*y+Cli(i));
            xtmpri_llane=(para(1,1)./(y-D_pre)+Bri(i)*y+Cri(i));
        else %model is line
            xtmpli_llane=(-eli(i)*y-fli(i)/2);
            xtmpri_llane=(-eri(i)*y-fri(i)/2);
        end
        xtmpli(i,:)=xtmpli_llane; xtmpri(i,:)=xtmpri_llane;
    end
    %convert to 480x640 coordinate
    xtmpli(i,:)=-(xtmpli(i,:)-cc_x_left);
    xtmpri(i,:)=-(xtmpri(i,:)-cc_x_right);
end
if leftrightlane==0 %left lane
    xtmpli_lb=xtmpli(1,:);
    xtmpli_rb=xtmpli(2,:);
    xtmpri_lb=xtmpri(1,:);
    xtmpri_rb=xtmpri(2,:);
else
    xtmpli_lb=xtmpli(2,:);
    xtmpli_rb=xtmpli(1,:);
    xtmpri_lb=xtmpri(2,:);
    xtmpri_rb=xtmpri(1,:);
end
xtmpl_imagine=0.5*(xtmpli_lb+xtmpli_rb);
xtmpr_imagine=0.5*(xtmpri_lb+xtmpri_rb);
ytmp_imagine=-(-y+cc_y+1)/2; %convert to -240x640

xtmpli_lb(xtmpli_lb<=1)=1;xtmpli_lb(xtmpli_lb>=ImgSize2)=ImgSize2;
xtmpli_rb(xtmpli_rb<=1)=1;xtmpli_rb(xtmpli_rb>=ImgSize2)=ImgSize2;
xtmpri_lb(xtmpri_lb<=1)=1;xtmpri_lb(xtmpri_lb>=ImgSize2)=ImgSize2;
xtmpri_rb(xtmpri_rb<=1)=1;xtmpri_rb(xtmpri_rb>=ImgSize2)=ImgSize2;

dumy2=C_in_region(dumy,xtmpli_lb,xtmpli_rb,xtmpri_lb+ImgSize2+30,xtmpri_rb+ImgSize2+30,ytopr,ybtmr);
dumy2(1:ytopr,:)=0; dumy2(ybtmr:end,:)=0;
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
    num_inlier=0;
else
    if sum(abs(para(:,2)))==0 %model is hypobola
%         [inlierl, inlierr, para2]=C_RANSAC_otherlane_hypo(leftrightlane, Datapointl, Datapointr, para(:,1), D_pre, Phai, H, cc_y);
        [inlierl, inlierr, para2]=C_RANSAC_otherlane_hypo_new(leftrightlane, Datapointl, Datapointr, para(:,1), D_pre, Phai, H, cc_y);
    else
%         [inlierl, inlierr, para2]=C_RANSAC_otherlane_line(leftrightlane, Datapointl, Datapointr, para(:,2), D_pre, Phai, H, cc_y);
        [inlierl, inlierr, para2]=C_RANSAC_otherlane_line_new(leftrightlane, Datapointl, Datapointr, para(:,2), D_pre, Phai, H, cc_y);
    end
    num_inlier=sum(inlierl)+sum(inlierr);
end


if num_inlier>50
    if sum(abs(para2(:,2)))==0 %model is hypobola
        xtmpl2=para2(1)./(y-D_pre)+y*para2(2)+para2(3);
        xtmpr2=para2(1)./(y-D_pre)+y*para2(4)+para2(5);
        lane_width=0.5*(abs((para2(2)-para(2))*H/cos(Phai))+abs((para2(3)-para(3))*H/(cos(Phai)*D_pre)));
    else %model is line
        xtmpl2=-para2(2,2)*y-0.5*para2(3,2);
        xtmpr2=-para2(4,2)*y-0.5*para2(5,2);
        lane_width=0.5*(abs((para2(2,2)-para(2,2))*H/cos(Phai))+abs((para2(3,2)-para(3,2))*H/(cos(Phai)*D_pre))/2.0);
    end
    xtmpl2=-xtmpl2+cc_x_left;
    ytmpl2=-(-y+cc_y+1)/2; %convert to -240x640
    xtmpr2=-xtmpr2+cc_x_right;
    ytmpr2=-(-y+cc_y+1)/2; %convert to -240x640
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
                sel_xtmpl2=[];sel_ytmpl2=[];sel_xtmpr2=[];sel_ytmpr2=[];lane_width=[];
                return;
            end
        end
    end
    
    ind1=intersectAB(-(-yl(inlierl==1)+cc_y+1)/2,ytmpl2);
    indytmp=find(ind1==1);
    tmp=floor(length(indytmp)/21);
    if tmp>=1
        equal_ind=(2:21)*tmp;
    else
        equal_ind=1:length(indytmp);
    end
    sel_xtmpl2=xtmpl2(indytmp(equal_ind));
    sel_ytmpl2=ytmpl2(indytmp(equal_ind));
    
    ind2=intersectAB(-(-yr(inlierr==1)+cc_y+1)/2,ytmpr2);
    indytmp=find(ind2==1);
    tmp=floor(length(indytmp)/21);
    if tmp>=1
        equal_ind=(2:21)*tmp;
    else
        equal_ind=1:length(indytmp);
    end
    sel_xtmpr2=xtmpr2(indytmp(equal_ind));
    sel_ytmpr2=ytmpr2(indytmp(equal_ind));  
end


