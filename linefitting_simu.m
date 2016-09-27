function [sampled_xtmpl, sampled_ytmpl, xtmpl, ytmpl, sampled_xtmpr, sampled_ytmpr, xtmpr, ytmpr, num_in, para]=linefitting_simu(kcsditmpl,kcsditmpr,th, D_pre, H_pre, num_ite)
% tic;
% cc_x_left=320; cc_x_right=349;
global_varibles;
ImgSize2=size(kcsditmpl,2);
[indvl, indul]=find(kcsditmpl>connected_pixel*th);
indvl=-(2*indvl-1)+cc_y;
indul=-indul+cc_x_left;
Datapointl=[indvl,indul]';
SubDatapointl=Datapointl;
[indvr, indur]=find(kcsditmpr>connected_pixel*th);
indvr=-(2*indvr-1)+cc_y;
indur=-indur+cc_x_right;
Datapointr=[indvr,indur]';
SubDatapointr=Datapointr;
%(indur,indvr) and (indul,indvl) all in camera coordinate (origin at
%principle point, x axis to left, y axis to up)
% toc;
% [line, bedf, tmp]=RANSAC(Datapoint, SubDatapoint);
% toc;
if (isempty(Datapointl) || isempty(Datapointr))
    linel=0; liner=0;
else
    [linel, liner, para]=P_Simu_RANSAC(Datapointl, SubDatapointl, Datapointr, SubDatapointr, D_pre, H_pre, num_ite, focal, dist_left_right, cc_y);
end
% sum(line)
% toc;
num_in=sum(linel)+sum(liner);
if sum(linel)==0 || sum(liner)==0 || sum(abs(para(:)))==0 || num_in<100
    xtmpl=[]; ytmpl=[]; sampled_xtmpl=[]; sampled_ytmpl=[];
    xtmpr=[]; ytmpr=[]; sampled_xtmpr=[]; sampled_ytmpr=[];
    para=[]; num_in=0;
    return;
end

top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even

if sum(abs(para(:,1)))==0
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
    equal_ind=(2:21)*tmp;
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
    equal_ind=(2:21)*tmp;
else
    equal_ind=1:length(indytmp);
end
sampled_xtmpr=xtmpr(indytmp(equal_ind));
sampled_ytmpr=ytmpr(indytmp(equal_ind));
end