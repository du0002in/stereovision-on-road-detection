function [sampled_xtmpl, sampled_ytmpl, xtmpl, ytmpl, sampled_xtmpr, sampled_ytmpr, xtmpr, ytmpr, num_in, para,inlier_xy_r, paraz, break_location,dumy_log]= ...
    linefitting_simu2(kcsditmpl,kcsditmpr, D_pre, H_pre, num_ite, pre_line_type, pre_k, zebra_status)
% tic;
% cc_x_left=320; cc_x_right=349;
global_varibles;
ImgSize2=size(kcsditmpl,2);
[indvl, indul]=find(kcsditmpl==1);
indvl=-(2*indvl-1)+cc_y;
indul=-indul+cc_x_left;
Datapointl=[indvl,indul]';
SubDatapointl=Datapointl;
[indvr, indur]=find(kcsditmpr==1);
indvr=-(2*indvr-1)+cc_y;
indur=-indur+cc_x_right;
Datapointr=[indvr,indur]';
SubDatapointr=Datapointr;
%(indur,indvr) and (indul,indvl) all in camera coordinate (origin at
%principle point, x axis to left, y axis to up)
% toc;
% [line, bedf, tmp]=RANSAC(Datapoint, SubDatapoint);
% toc;
% if (isempty(Datapointl) || isempty(Datapointr))
if length(Datapointl)<4 || length(Datapointr)<4
    linel=0; liner=0;linelz=0;linerz=0;para=zeros(5,2);paraz=zeros(5,2);
else
    if zebra_status==0
        [linel, liner, para]=P_Simu_RANSAC(Datapointl, SubDatapointl, Datapointr, SubDatapointr, D_pre, H_pre, num_ite, focal, dist_left_right, cc_y, pre_k);
    else
        [linel, liner, para]=P_Zebra_RANSAC_single_line(Datapointl, SubDatapointl, Datapointr, SubDatapointr, D_pre, H_pre, num_ite, focal, dist_left_right, cc_y);
        para=[zeros(5,1) para];
    end
    [linelz, linerz, paraz, y_int]=P_Zebra_RANSAC_new(Datapointl, SubDatapointl, Datapointr, SubDatapointr, D_pre, H_pre, num_ite, focal, dist_left_right, cc_y);
end
% sum(line)
% toc;
C0=4*para(1,1)*cos(atan(D_pre/focal))^3/H_pre/focal^2;
C0_z=0;
sig_prob=0.04;
prob=1/2.506628/sig_prob*exp(-((C0-pre_k)-1.4854e-5)^2/2/sig_prob^2);
if prob<5
    prob=5;
end
prob_z=1/2.506628/sig_prob*exp(-((C0_z-pre_k)-1.4854e-5)^2/2/sig_prob^2);
if prob_z<5
    prob_z=5;
end

num_in=sum(linel)+sum(liner);
num_inz=sum(linelz)+sum(linerz);
dumy_log=0;
% a=sum(linel); b=sum(liner);
if (sum(linel)<30 || sum(liner)<30 || sum(abs(para(:)))==0 || num_in<60) && ...         %both cannot
        (sum(linelz)<30 || sum(linerz)<30 || sum(abs(paraz(:)))==0 || num_inz<60)       %both cannot
    xtmpl=[]; ytmpl=[]; sampled_xtmpl=[]; sampled_ytmpl=[];
    xtmpr=[]; ytmpr=[]; sampled_xtmpr=[]; sampled_ytmpr=[];
    para=zeros(5,2); num_in=0; inlier_xy_r=[]; break_location=[];paraz=zeros(5,2);
    return;
elseif (sum(linel)<30 || sum(liner)<30 || sum(abs(para(:)))==0 || num_in<60) && ...     %normal line cannot
        ~(sum(linelz)<30 || sum(linerz)<30 || sum(abs(paraz(:)))==0 || num_inz<60)      %zebra can
    if (sum(indvl(linelz==1)>y_int(1))+sum(indvr(linerz==1)>y_int(2)))>30 && ...
            (sum(indvl(linelz==1)<=y_int(1))+sum(indvr(linerz==1)<=y_int(2)))>30
        para=paraz;
        linel=linelz;
        liner=linerz;
        num_in=num_inz;
        y_int_ave=mean(y_int);
%         dumy_log=[sum(indvl(linelz==1)>y_int(1))+sum(indvr(linerz==1)>y_int(2)) sum(indvl(linelz==1)<=y_int(1))+sum(indvr(linerz==1)<=y_int(2))];
    else
        xtmpl=[]; ytmpl=[]; sampled_xtmpl=[]; sampled_ytmpl=[];
        xtmpr=[]; ytmpr=[]; sampled_xtmpr=[]; sampled_ytmpr=[];
        para=zeros(5,2); num_in=0; inlier_xy_r=[]; break_location=[];paraz=zeros(5,2);
        return;
    end
elseif ~(sum(linel)<30 || sum(liner)<30 || sum(abs(para(:)))==0 || num_in<60) && ...    %normal line can
        (sum(linelz)<30 || sum(linerz)<30 || sum(abs(paraz(:)))==0 || num_inz<60)       %zebra cannot
    paraz=zeros(5,2);break_location=[];
else                                                                                    %both can
    if num_inz*prob_z<=0.8*num_in*prob
        paraz=zeros(5,2);break_location=[];
    elseif num_inz*prob_z>1.3*num_in*prob && (sum(indvl(linelz==1)>y_int(1))+sum(indvr(linerz==1)>y_int(2)))>30 && ...
            (sum(indvl(linelz==1)<=y_int(1))+sum(indvr(linerz==1)<=y_int(2)))>30
        para=paraz;
        linel=linelz;
        liner=linerz;
        num_in=num_inz;
        y_int_ave=mean(y_int);
%         dumy_log=[sum(indvl(linelz==1)>y_int(1))+sum(indvr(linerz==1)>y_int(2)) sum(indvl(linelz==1)<=y_int(1))+sum(indvr(linerz==1)<=y_int(2))];
    else
        if (pre_line_type==3) && (sum(indvl(linelz==1)>y_int(1))+sum(indvr(linerz==1)>y_int(2)))>30 && ...
            (sum(indvl(linelz==1)<=y_int(1))+sum(indvr(linerz==1)<=y_int(2)))>30
            para=paraz;
            linel=linelz;
            liner=linerz;
            num_in=num_inz;
            y_int_ave=mean(y_int);
%             dumy_log=[sum(indvl(linelz==1)>y_int(1))+sum(indvr(linerz==1)>y_int(2)) sum(indvl(linelz==1)<=y_int(1))+sum(indvr(linerz==1)<=y_int(2))];
        else
            paraz=zeros(5,2);break_location=[];
        end
    end
end



top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even

if sum(paraz(:))==0
    if (sum(abs(para(:,1)))==0)
        ytmpl=top_row:-2:btm_row;
        xtmpl=(-para(3,2)-2*para(2,2)*ytmpl)/(2*para(1,2));
        ytmpr=top_row:-2:btm_row;
        xtmpr=(-para(5,2)-2*para(4,2)*ytmpr)/(2*para(1,2));
    else
        if mod(floor(D_pre),2)==0
            D_pre=floor(D_pre)-1;
        end
        ytmpl=min(floor(D_pre),top_row):-2:btm_row;
        xtmpl=para(1,1)./(ytmpl-D_pre)+para(2,1)*ytmpl+para(3,1);
        ytmpr=min(floor(D_pre),top_row):-2:btm_row;
        xtmpr=para(1,1)./(ytmpr-D_pre)+para(4,1)*ytmpr+para(5,1);
    end
else
    ytmpl=top_row:-2:btm_row;
    xtmpl=(-para(3,2)-2*para(2,2)*ytmpl(ytmpl>=y_int_ave))/(2*para(1,2));
    xtmpl=[xtmpl (-para(3,1)-2*para(2,1)*ytmpl(ytmpl<y_int_ave))/(2*para(1,1))];
    ytmpr=top_row:-2:btm_row;
    xtmpr=(-para(5,2)-2*para(4,2)*ytmpr(ytmpr>=y_int_ave))/(2*para(1,2));
    xtmpr=[xtmpr (-para(5,1)-2*para(4,1)*ytmpr(ytmpr<y_int_ave))/(2*para(1,1))];
    break_location=find(ytmpr<y_int_ave,1);
    y_int_ave=-(-y_int_ave+cc_y+1)/2;
end

dumy_log=.413*684/(xtmpr(end)-xtmpl(end));

% toc
ytmpl=-(-ytmpl+cc_y+1)/2;
xtmpl=-(xtmpl-cc_x_left);
ytmpr=-(-ytmpr+cc_y+1)/2;
xtmpr=-(xtmpr-cc_x_right);

ytmpl(round(xtmpl)<=0)=[];
xtmpl(round(xtmpl)<=0)=[];
ytmpl(round(xtmpl)>ImgSize2)=[];
xtmpl(round(xtmpl)>ImgSize2)=[];

ytmpr(round(xtmpr)<=0)=[];
xtmpr(round(xtmpr)<=0)=[];
ytmpr(round(xtmpr)>ImgSize2)=[];
xtmpr(round(xtmpr)>ImgSize2)=[];

ind1=intersectAB(-(-indvl(linel==1)+cc_y+1)/2,ytmpl);
% ind1=intersectAB(ytmpl,ytmpl);
indytmp=find(ind1==1);
tmp=floor(length(indytmp)/21);
if tmp>=1
    equal_ind=round((2:21)*length(indytmp)/21);
else
    equal_ind=1:length(indytmp);
end
sampled_xtmpl=xtmpl(indytmp(equal_ind));
sampled_ytmpl=ytmpl(indytmp(equal_ind));

ind2=intersectAB(-(-indvr(liner==1)+cc_y+1)/2,ytmpr);
% ind2=intersectAB(ytmpr,ytmpr);
indytmp=find(ind2==1);
tmp=floor(length(indytmp)/21);
if tmp>=1
    equal_ind=round((2:21)*length(indytmp)/21);
else
    equal_ind=1:length(indytmp);
end
sampled_xtmpr=xtmpr(indytmp(equal_ind));
sampled_ytmpr=ytmpr(indytmp(equal_ind));
inlier_xy_r=[-(indur(liner==1)-cc_x_right) -(-indvr(liner==1)+cc_y+1)/2];
if sum(paraz(:))~=0
    break_location=find(ytmpr<y_int_ave,1);
end

end