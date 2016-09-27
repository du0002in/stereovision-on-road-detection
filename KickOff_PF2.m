ii=ite_s;
rr_arrow=[]; rr_tip=[];log_p_tip(ii)=0;log_zebra(ii)=0;log_hump(ii)=0;
num_in_pre=0; para=zeros(3,2);
% Initilize particles
Particle=zeros(8,P_No);
Particle(3,:)=pi/2; %theta is 90 degree
% [v_m, phai_m]=get_measurements;
v_m=act_v(ii); phai_m=act_phai(ii);
if v_m~=0
    v0=v_m+sqrt(0.1).*randn(1,P_No);
else
    v0=zeros(1,P_No);
end
phai0=phai_m+sqrt(0.1).*randn(1,P_No);
Particle(4,:)=cos(Particle(3,:)).*v0;
Particle(5,:)=sin(Particle(3,:)).*v0;
Particle(6,:)=tan(phai0).*v0/L_wheels;
Particle(7,:)=v0;
Particle(8,:)=phai0;
log_hump(ii)=0;
log_zebra(ii)=0;
in_start=tic;
leftimg=imread(strcat(leftimg_dir,int2str(ii),'.jpg'));
rightimg=imread(strcat(rightimg_dir,int2str(ii),'.jpg'));
start=tic;
leftimg=double(max(leftimg,[],3));
rightimg=double(max(rightimg,[],3));
% tic
leftimg=Retify(leftimg,a1_left,a2_left,a3_left,a4_left,ind_1_left,ind_2_left,ind_3_left,ind_4_left,ind_new_left);
rightimg=Retify(rightimg,a1_right,a2_right,a3_right,a4_right,ind_1_right,ind_2_right,ind_3_right,ind_4_right,ind_new_right);
img=[leftimg(1:2:ImgSize1,:) zeros(ImgSize1/2, 30) rightimg(1:2:ImgSize1,:)];
img(235:end,:)=0;
mask=img;mask(img<230)=0;mask(img>=230)=1;
for ini_ridge_loop=1:5
    Ridge_width_ini=Ridge_width_ini+2.0;
    RidgeParameter;
    ImgSmooth=zeros(size(img));
    for i=1:(length(ix)-1)
        ImgSmooth(ix(i):(ix(i+1)-1),:)=conv2(img(ix(i):(ix(i+1)-1),:),xfil2{i},'same');
    end
    
    maxconvimg=C_new_Ridge(ImgSmooth,tdx,Offset);
    C_refine_Ridge(maxconvimg,Offset);
%     maxconvimg=maxconvimg.*mask;
%     ImgSmooth=conv2(ImgSmooth,yfil2,'same');
%     [kcsditmp,th]=C_Ridge(ImgSmooth,306400);
%     for i=1:8
%         convimg(:,:,i)=conv2(kcsditmp,wd(:,:,i),'same');
%     end
%     maxconvimg=max(convimg,[],3);
    dumy1tmp=maxconvimg;
%     dumy1tmp(dumy1tmp<connected_pixel*th)=0; dumy1tmp(dumy1tmp>=connected_pixel*th)=1;
    num_ite=35000; %number of iterations for RANSAC
    intensity_td0=70; %for night vision
%     intensity_td0=230; %for daytime vision
    for intensity_td_loop=1:5
        intensity_td0=intensity_td0+10;
        mask=img;mask(img<intensity_td0)=0;mask(img>=intensity_td0)=1;
        %          mask(:,[1:186, 670:944])=0;
        maxconvimg=maxconvimg.*mask;
        for pre_C0=[-0.05 0 0.05]
            [sel_xtmpl0,sel_ytmpl0, xtmpl0, ytmpl0,sel_xtmpr0,sel_ytmpr0, xtmpr0, ytmpr0, num_in0, para0, inlier_xy_r0,paraz0, break_location0]= ...
                linefitting_simu2(maxconvimg(:,1:ImgSize2),maxconvimg(:,end-ImgSize2+1:end),D_pre, H, num_ite,pre_line_type,pre_C0,0);
            if (sum(abs(para(:,1)))==0 && sum(abs(para0(:,1)))==0 && num_in0>num_in_pre) || ...
                    (sum(abs(para(:,2)))==0 && sum(abs(para0(:,2)))==0 && num_in0>num_in_pre) || ...
                    (sum(abs(para(:,1)))==0 && sum(abs(para0(:,2)))==0 && 0.8*num_in0>num_in_pre) || ...
                    (sum(abs(para(:,2)))==0 && sum(abs(para0(:,1)))==0 && num_in0>0.8*num_in_pre)
                dumy1=dumy1tmp;
                dumy2=maxconvimg;
                %         dumy2(dumy2<connected_pixel*th)=0;dumy2(dumy2>=connected_pixel*th)=1;
                dumy=[dumy1;dumy2];
                imwrite(dumy,strcat(save_dir,'\Ridge',int2str(ii),'.jpg'),'jpg');
                sel_xtmpl=sel_xtmpl0; sel_ytmpl=sel_ytmpl0; xtmpl=xtmpl0; ytmpl=ytmpl0;
                sel_xtmpr=sel_xtmpr0; sel_ytmpr=sel_ytmpr0; xtmpr=xtmpr0; ytmpr=ytmpr0;
                num_in=num_in0; para=para0; inlier_xy_r=inlier_xy_r0;
                num_in_pre=num_in0;
                paraz=paraz0; break_location=break_location0;
%                 log_tdx{ii}=tdx;
            end
        end
    end
end

% lognum_in(ii)=num_in;
if num_in<60
    diff_marking_width=0;
    logPhai(ii)=0.2757;
    logH(ii)=1.665;
    log_marking_width(ii)=marking_width_pre;
    log_marking_width_f(ii)=marking_width_pre;
    % % %         log_marking_width2(ii)=marking_width_pre2;
    endt(ii)=toc(start);
    endf(ii)=toc(in_start);
    disp 'Cannot initialize';
    return;
end

sel_xtmpl_tmp=round(sel_xtmpl); sel_xtmpl_tmp(sel_xtmpl_tmp<1)=1; sel_xtmpl_tmp(sel_xtmpl_tmp>ImgSize2)=ImgSize2;
sel_xtmpr_tmp=round(sel_xtmpr); sel_xtmpr_tmp(sel_xtmpr_tmp<1)=1; sel_xtmpr_tmp(sel_xtmpr_tmp>ImgSize2)=ImgSize2;
indtmp_l=sub2ind(size(img),-sel_ytmpl,sel_xtmpl_tmp);
intensity_l(ii,:)=[min(img(indtmp_l)) max(img(indtmp_l))];
indtmp_r=sub2ind(size(img),-sel_ytmpr,sel_xtmpr_tmp+ImgSize2+30);
intensity_r(ii,:)=[min(img(indtmp_r)) max(img(indtmp_r))];
intesnity_td_min=min([min(img(indtmp_l)) min(img(indtmp_r))])-5;
intesnity_td_max=max([max(img(indtmp_l)) max(img(indtmp_r))])+5;
%convert back to original image coordinate 480x640
% % [~,ir,il]=intersect(sel_ytmpr, ytmpl, 'stable');
if xtmpr(end)<320 %This condition may not be true...Need to update this based on Xworld coordinate
    leftrightlane=0;
else
    leftrightlane=1;
end
% % % % % if leftrightlane==0 %left lane
% % % % %     ir=intersectAB(ytmpl,sel_ytmpr);
% % % % %     il=intersectAB(sel_ytmpr,ytmpl);
% % % % %     % toc(start)
% % % % %     sel_ytmpr=sel_ytmpr(ir==1);
% % % % %     sel_xtmpr=sel_xtmpr(ir==1);
% % % % %     Corespx=xtmpl(il==1);
% % % % %     delta_x=round(0.308*(-sel_ytmpr)+7.385);
% % % % %     sel_ytmpr=2*sel_ytmpr-1; %convert back to -480x640 coordinate
% % % % %     ri_lb=round(sel_xtmpr-delta_x); %right image left boundary
% % % % %     ri_rb=round(sel_xtmpr+delta_x); %right image right boundary
% % % % %     ri_lb(ri_lb<=1)=1;
% % % % %     ri_rb(ri_rb>=ImgSize2)=ImgSize2;
% % % % %     li_lb=round(Corespx-2*delta_x);
% % % % %     li_rb=round(Corespx+2*delta_x);
% % % % %     li_lb(li_lb<=1)=1;
% % % % %     li_rb(li_rb>=ImgSize2)=ImgSize2;
% % % % %     ri_lb=ri_lb-1;
% % % % %     ri_rb=ri_rb-1;
% % % % %     li_lb=li_lb-1;
% % % % %     li_rb=li_rb-1;
% % % % %     
% % % % %     IntensityR=rightimg(-sel_ytmpr,:);
% % % % %     CorespL=leftimg(-sel_ytmpr,:);
% % % % %     % toc(start)
% % % % %     disparitytmp=ZNCC(ri_lb', ri_rb', li_lb', li_rb', IntensityR, CorespL);
% % % % %     % toc(start)
% % % % %     li_lb=li_lb+1;
% % % % %     ri_lb=ri_lb+1;
% % % % %     % Refer to the sign convention in 3D Vision notes
% % % % %     % Convert to left and right camera coordinates respectively
% % % % %     CorepxL=-(li_lb+disparitytmp'+(round(sel_xtmpr)-ri_lb)-cc_x_left);
% % % % %     CorepxR=-(round(sel_xtmpr)-cc_x_right);
% % % % % else %right lane
% % % % %     il=intersectAB(ytmpr,sel_ytmpl);
% % % % %     ir=intersectAB(sel_ytmpl,ytmpr);
% % % % %     % toc(start)
% % % % %     sel_ytmpl=sel_ytmpl(il==1);
% % % % %     sel_xtmpl=sel_xtmpl(il==1);
% % % % %     Corespx=xtmpr(ir==1);
% % % % %     delta_x=round(0.308*(-sel_ytmpl)+7.385);
% % % % %     sel_ytmpl=2*sel_ytmpl-1; %convert back to -480x640 coordinate
% % % % %     li_lb=round(sel_xtmpl-delta_x); %right image left boundary
% % % % %     li_rb=round(sel_xtmpl+delta_x); %right image right boundary
% % % % %     li_lb(li_lb<=1)=1;
% % % % %     li_rb(li_rb>=ImgSize2)=ImgSize2;
% % % % %     ri_lb=round(Corespx-2*delta_x);
% % % % %     ri_rb=round(Corespx+2*delta_x);
% % % % %     ri_lb(ri_lb<=1)=1;
% % % % %     ri_rb(ri_rb>=ImgSize2)=ImgSize2;
% % % % %     li_lb=li_lb-1;
% % % % %     li_rb=li_rb-1;
% % % % %     ri_lb=ri_lb-1;
% % % % %     ri_rb=ri_rb-1;
% % % % %     
% % % % %     IntensityL=leftimg(-sel_ytmpl,:);
% % % % %     CorespR=rightimg(-sel_ytmpl,:);
% % % % %     % toc(start)
% % % % %     disparitytmp=ZNCC(li_lb', li_rb', ri_lb', ri_rb', IntensityL, CorespR);
% % % % %     % toc(start)
% % % % %     li_lb=li_lb+1;
% % % % %     ri_lb=ri_lb+1;
% % % % %     % Refer to the sign convention in 3D Vision notes
% % % % %     % Convert to left and right camera coordinates respectively
% % % % %     CorepxR=-(ri_lb+disparitytmp'+(round(sel_xtmpl)-li_lb)-cc_x_right);
% % % % %     CorepxL=-(round(sel_xtmpl)-cc_x_left);
% % % % % end
if leftrightlane==0 %left lane
    [~,il,ir]=intersect(ytmpl,sel_ytmpr,'stable');
    sel_ytmpr=sel_ytmpr(ir);
    sel_xtmpr=sel_xtmpr(ir);
    CorepxR=-(sel_xtmpr-cc_x_right);
    CorepxL=-(xtmpl(il)-cc_x_left);
    sel_ytmpr=2*sel_ytmpr-1; %convert back to -480x640 coordinate
else
    [~,ir,il]=intersect(ytmpr,sel_ytmpl,'stable');
    sel_ytmpl=sel_ytmpl(il);
    sel_xtmpl=sel_xtmpl(il);
    CorepxL=-(sel_xtmpl-cc_x_left);
    CorepxR=-(xtmpr(ir)-cc_x_right);
    sel_ytmpl=2*sel_ytmpl-1; %convert back to -480x640 coordinate
end
disparity=CorepxR-CorepxL;
indtmp=(disparity~=0);
Zc=dist_left_right*focal./disparity;
Xc=CorepxL.*Zc/focal+dist_left_right/2;
% or Xc=CorepxR.*Zc/focal-dist_left_right/2;
if leftrightlane==0
    Yc=(round(sel_ytmpr)+cc_y).*Zc/focal;
    ind_near=(sel_ytmpr<-192); %Use those pixels at 2/3bottom to calculate phai
else
    Yc=(round(sel_ytmpl)+cc_y).*Zc/focal;
    ind_near=(sel_ytmpl<-192); %Use those pixels at 2/3bottom to calculate phai
end
ind_near=ind_near & indtmp;
Xc1=Xc; Zc1=Zc; Yc1=Yc;
marking_length(ii)=((Xc(1)-Xc(end))^2+(Zc(1)-Zc(end))^2)^0.5;
%     if length(Yc)>5
%         Yc_near=Yc(end-5:end); Zc_near=Zc(end-5:end);
%     else
%         Yc_near=Yc; Zc_near=Zc;
%     end
if sum(ind_near)>5
    Yc_near=Yc(ind_near); Zc_near=Zc(ind_near); Xc_near=Xc(ind_near);
else
    Yc_near=Yc; Zc_near=Zc; Xc_near=Xc;
end
Atmp=[Zc_near' -ones(length(Zc_near),1)];
Btmp=Yc_near'*cos(Theta)+Xc_near'*sin(Theta);
tan_Phai_H=pinv(Atmp)*Btmp;
tan_Phai=tan_Phai_H(1);
Phai=atan(tan_Phai);
logPhaiNaN(ii)=Phai;
if isnan(Phai)
    Phai=atan(D_pre/focal);
end
logPhai(ii)=Phai;
Phai=median(logPhai(max([ite_s, ii-10]):ii)); %median value of the previous 10 results.
logPhai_f(ii)=Phai;
D_pre=focal*tan(Phai);
%     H=mean((Zc_near*tan(Phai)-Yc_near)*cos(Phai));
H=tan_Phai_H(2)*cos(Phai);
logHNaN(ii)=H;
if H<1.4 || isnan(H) || H>1.8
    H=1.665;
end
logH(ii)=H;
H=median(logH(max([ite_s, ii-10]):ii)); %median value of the previous 10 results.
logH_f(ii)=H;
% [xtmpl, ytmpl, xtmpr, ytmpr, refine_para]=refine_model(leftrightlane,para,sel_ytmpr,sel_ytmpl,D_pre,CorepxR,CorepxL);
refine_para=para;
if numel(ytmpl)==0 || numel(ytmpr)==0
    endt(ii)=toc(start);
    endf(ii)=toc(in_start);
    disp 'Cannot initialize';
    return;
end
if sum(paraz(:)==0)
    [sel_xtmpl2,sel_ytmpl2,xtmpl2, ytmpl2, sel_xtmpr2,sel_ytmpr2,xtmpr2, ytmpr2, num_in2, para2,lwtmp,~,xtmpl_imagine, xtmpr_imagine, ytmp_imagine]= ...
        fit_other_lane(leftrightlane,dumy1,refine_para,-ytmpr(1),-ytmpr(end),H,Phai,D_pre, ximg, yimg);
    if leftrightlane==0
        expos_box={xtmpl, xtmpl_imagine, xtmpr, xtmpr_imagine, -ytmpr(1), 240, ytmpr, ytmp_imagine};
    else
        expos_box={xtmpl_imagine, xtmpl, xtmpr_imagine, xtmpr, -ytmpr(1), 240, ytmp_imagine, ytmpr};
    end
else
    xtmpl2=[];ytmpl2=[];xtmpr2=[];ytmpr2=[];
    sel_xtmpl2=[];sel_ytmpl2=[];sel_xtmpr2=[];sel_ytmpr2=[];
    num_in2=0; para2=[]; lwtmp=[];
    if leftrightlane==0
        expos_box={xtmpl, xtmpl+round(46.8-2.26*ytmpl), xtmpr, xtmpr+round(46.8-2.26*ytmpr), -ytmpr(1), 240, ytmpr, ytmpr};
    else
        expos_box={xtmpl-round(46.8-2.26*ytmpl), xtmpl, xtmpr-round(46.8-2.26*ytmpr), xtmpr, -ytmpr(1), 240, ytmpr, ytmpr};
    end
end

if isempty(lwtmp)
    log_lw(ii,1)=0;
else
    log_lw(ii,1)=lwtmp;
end
log_lane_width_f(ii)=lane_width;
R=[cos(Theta) -sin(Theta) 0; ...
    cos(Phai)*sin(Theta) cos(Phai)*cos(Theta) -sin(Phai); ...
    sin(Phai)*sin(Theta) sin(Phai)*cos(Theta) cos(Phai);];
T=[0;H;L_wheel_cam];
World=[R T+offset_vector;0 0 0 1]*[Xc; Yc; Zc; ones(1,length(Zc))];
% logXw=World(1); logYw=World(2); logZw=World(3);
% toc
Xw=World(1,:); Yw=World(2,:); Zw=World(3,:);
inf_ind=(isinf(Xw) | isinf(Yw) | isinf(Zw))| (isnan(Xw) | isnan(Yw) | isnan(Zw));
Xw(inf_ind)=[];Yw(inf_ind)=[];Zw(inf_ind)=[];
World(:,inf_ind)=[];

World1=[R T+offset_vector;0 0 0 1]*[Xc1; Yc1; Zc1; ones(1,length(Zc1))];
Xw1=World1(1,:); Yw1=World1(2,:); Zw1=World1(3,:);
inf_ind1=(isinf(Xw1) | isinf(Yw1) | isinf(Zw1)) | (isnan(Xw1) | isnan(Yw1) | isnan(Zw1)) | (Zw1>20);
Xw1(inf_ind1)=[];Yw1(inf_ind1)=[];Zw1(inf_ind1)=[];
World1(:,inf_ind1)=[];
World2=World1;
World3=World2;

if sum(paraz(:)==0)
    if sum(abs(para(:,2)))==0 %model is hypobola
        Atmp=[Zw.^2; Zw; ones(1,length(Zw))]';
        Btmp=Xw';
        abc=pinv(Atmp)*Btmp;
        pre_line_type=2;
    else %model is line
        Atmp=[Zw; ones(1,length(Zw))]';
        Btmp=Xw';
        bc=pinv(Atmp)*Btmp;
        abc=[0; bc];
        pre_line_type=1;
    end
else
    [Xw, Yw, Zw, abc]=equavalent_line(break_location, Xw, Zw);
    World=[Xw;Yw;Zw];
    log_line_para(ii)=3;
    pre_line_type=3;
end
[xc_m, d_theta_m]=solve_xc_d_theta(abc,0,0,pi/2);
pre_abc1=abc;
pre_World=World;
pre_para=para;
% log_Yw_mean(ii)=mean(Yw);

log_xc_m(ii)=xc_m;
log_d_theta_m(ii)=d_theta_m;
if num_in2>num_in
%     xtmpl=xtmpl2; ytmpl=ytmpl2; xtmpr=xtmpr2; ytmpr=ytmpr2;
    lanelinechange=1;
else lanelinechange=0;
end
%     logWidth(ii)=f_width;
lognum_in2(ii)=num_in2;
if sum(ind_near)>6
    ind_near(1:(end-7))=0;
end
if mod(sum(ind_near),2)==0
    indtmp=find(ind_near==1,1);
    ind_near(indtmp)=0;
end
if sum(ind_near)>=5
    if leftrightlane==0
        y_for_tdx=sel_ytmpr(ind_near)';
        x_for_tdx=round(sel_xtmpr(ind_near)');
        Intensity=rightimg(-y_for_tdx,:);
    else
        y_for_tdx=sel_ytmpl(ind_near)';
        x_for_tdx=round(sel_xtmpl(ind_near)');
        Intensity=leftimg(-y_for_tdx,:);
    end
    y_for_tdx=y_for_tdx+cc_y; %convert to camera coordinate
    tdxtmp1=round((cos(Phai)/H*(D_pre-y_for_tdx)*0.10)/2);
    tdxtmp2=round((cos(Phai)/H*(D_pre-y_for_tdx)*0.15)/2);
    tdxtmp3=round((cos(Phai)/H*(D_pre-y_for_tdx)*0.20)/2);
    quater_tdxtmp3=round(1/4*tdxtmp3);
    x_for_tdx_lb=x_for_tdx-tdxtmp3-quater_tdxtmp3;
    x_for_tdx_rb=x_for_tdx+tdxtmp3+quater_tdxtmp3;
    x_for_tdx_lb(x_for_tdx_lb<1)=1;
    x_for_tdx_lb(x_for_tdx_lb>ImgSize2)=ImgSize2;
    x_for_tdx_rb(x_for_tdx_rb<1)=1;
    x_for_tdx_rb(x_for_tdx_rb>ImgSize2)=ImgSize2;
    marking_width_all=C_find_lane_marking_width(Intensity,x_for_tdx_lb,x_for_tdx_rb,tdxtmp1,tdxtmp2,tdxtmp3,x_for_tdx);
    marking_width=median(marking_width_all);
    
    log_marking_width_left(ii)=marking_width;
    marking_width=median(log_marking_width_left(max([ite_s, ii-4]):ii));
    log_marking_width_f_left(ii)=marking_width;
    log_marking_width_right(ii)=marking_width;
    log_marking_width_f_right(ii)=marking_width;
    
%     if marking_width~=marking_width_pre
        marking_width_pre=marking_width;
        if mod(floor(D_pre),2)==0
            D_pre1=floor(D_pre)-1;
        else D_pre1=floor(D_pre);
        end
        yy=min(floor(D_pre1),top_row):-1:btm_row;
        tdxtmp=(cos(Phai)/H*(D_pre-yy)*marking_width)/2; %0.10m is the lane marking width
        if length(tdxtmp)>=length(tdx)
            tdx=round(tdxtmp((end-length(tdx)+1):end));
        else
            tdx=round([0.5*ones(1,(length(tdx)-length(tdxtmp))) tdxtmp]);
        end
        tdx(tdx==0)=1;
        [~,ix,~]=unique(tdx,'stable');ix=ix';
        j=1;
        for i=ix
            weight3=exp((-0.5*(-tdx(i):tdx(i)).^2/(tdx(i))^2));
            weight3=weight3/sum(weight3);
            xfil2{j}=weight3;
            j=j+1;
        end
        ix(length(ix)+1)=ImgSize1-Offset+1;
        ix=round(ix/2)+Offset/2;
%     end
else
    log_marking_width_left(ii)=marking_width_pre;
    log_marking_width_f_left(ii)=marking_width_pre;
    log_marking_width_right(ii)=marking_width_pre;
    log_marking_width_f_right(ii)=marking_width_pre;
end
endt(ii)=toc(start);
endf(ii)=toc(in_start);


log_Estimate_xc(ii)=xc_m;
log_Estimate_d_theta(ii)=d_theta_m;
Estimate_d_theta=d_theta_m;

if leftrightlane==0
    log_xc_m(ii)=log_xc_m(ii)-lane_width;
    log_Estimate_xc(ii)=log_Estimate_xc(ii)-lane_width;
%     log_predi_xc(ii)=log_predi_xc(ii)-lane_width;
elseif leftrightlane==2
    log_xc_m(ii)=log_xc_m(ii)+lane_width;
    log_Estimate_xc(ii)=log_Estimate_xc(ii)+lane_width;
%     log_predi_xc(ii)=log_predi_xc(ii)+lane_width;
end
logleftrightlane(ii)=leftrightlane;

dt=act_t(ii+1);
dt=0.2;
Predict_P=zeros(8,P_No);
Predict_P(4,:)=cos(Particle(3,:)).*Particle(7,:);
Predict_P(5,:)=sin(Particle(3,:)).*Particle(7,:);
Predict_P(6,:)=-tan(Particle(8,:)).*Particle(7,:)/L_wheels;
Predict_P(1,:)=Particle(1,:)+Particle(4,:)*dt+0.25*dt/0.2*marking_width*randn(1,P_No);
Predict_P(2,:)=Particle(2,:)+Particle(5,:)*dt+0.25*dt/0.2*marking_width*randn(1,P_No);
Predict_P(3,:)=Particle(3,:)+Particle(6,:)*dt;

Predict_P(7,:)=Particle(7,:);
Predict_P(8,:)=Particle(8,:);
[Predict_Xc, Predict_d_theta]=C_solve_xc_d_theta(abc, Predict_P(1:3,:));

loop_t=logloop_data{ii}{1};
loop_phi=logloop_data{ii}{2};
loop_v=logloop_data{ii}{3};
RT_abc=GetRotMat([0;0;pi/2], loop_t, loop_phi, loop_v);
ZZ=0:0.2:4; XX=pre_abc1(1)*ZZ.^2+pre_abc1(2)*ZZ+pre_abc1(3);
rot_World=RT_abc*[XX;ZZ;ones(1,length(ZZ))];

particle_index=1:P_No;
% 
% if xc_m>=0
%     leftrightlane=0;
% else
%     leftrightlane=1;
% end
pre_leftrightlane=leftrightlane;
if lanelinechange==1
    leftrightlane=abs(leftrightlane-1);
end

if sum(paraz(:))==0
    logC0(ii)=4*para(1,1)*cos(Phai)^3/H/focal^2;
else logC0(ii)=0;
end

plotimg;
% % % if leftrightlane==0
% % %     log_xc_m(ii)=log_xc_m(ii)-lane_width;
% % %     log_Estimate_xc(ii)=log_Estimate_xc(ii)-lane_width;
% % % %     log_predi_xc(ii)=log_predi_xc(ii)-lane_width;
% % % end
% % % logleftrightlane(ii)=leftrightlane;

% % % % % Predict_x_mean=mean(Predict_P(1,:));
% % % % % Predict_z_mean=mean(Predict_P(2,:));
% % % % % Predict_theta_mean=mean(Predict_P(3,:));
% % % % % Rtmp2=[cos(Predict_theta_mean-0.5*pi) sin(Predict_theta_mean-0.5*pi); ...
% % % % %         -sin(Predict_theta_mean-0.5*pi) cos(Predict_theta_mean-0.5*pi)];
% newWorld=[cos(Predict_theta_mean-0.5*pi) sin(Predict_theta_mean-0.5*pi) -Predict_x_mean; ...
%         -sin(Predict_theta_mean-0.5*pi) cos(Predict_theta_mean-0.5*pi) -Predict_z_mean;]* ...
%         World([1,3,4],:);
% % % % % newWorld=[Rtmp2 -Rtmp2*[Predict_x_mean;Predict_z_mean]]*World([1,3,4],:);
% % % % % newWorld=[newWorld(1,:); World(2,:); newWorld(2,:); World(4,:)];
% % % % % newXYZc=([R T;0 0 0 1])\newWorld;
% % % % % newXc=newXYZc(1,:); newYc=newXYZc(2,:); newZc=newXYZc(3,:);
% % % % % yimg=focal*newYc./newZc; %[cc_y (cc_y-480)]*640
% % % % % ximgl=(newXc-dist_left_right/2)*focal./newZc;
% % % % % ximgr=(newXc+dist_left_right/2)*focal./newZc;
% % % % % if sum(abs(para(:,1)))==0
% % % % %     Atmp=-[yimg' ones(length(yimg),1)]; Btmp=ximgl';
% % % % %     ef_l=pinv(Atmp)*Btmp;
% % % % %     Btmp=ximgr';
% % % % %     ef_r=pinv(Atmp)*Btmp;
% % % % %     ytmpl=top_row:-2:btm_row;
% % % % %     xtmpl=-(-ef_l(1)*ytmpl-ef_l(2))+cc_x_left;
% % % % %     ytmpr=top_row:-2:btm_row;
% % % % %     xtmpr=-(-ef_r(1)*ytmpr-ef_r(2))+cc_x_right;
% % % % % else
% % % % %     Atmp=[1./(yimg'-D_pre) yimg' ones(length(yimg),1)]; Btmp=ximgl';
% % % % %     ABC_l=pinv(Atmp)*Btmp;
% % % % %     Btmp=ximgr';
% % % % %     ABC_r=pinv(Atmp)*Btmp;
% % % % %     ytmpl=min(floor(D_pre),top_row):-2:btm_row;
% % % % %     xtmpl=-(ABC_l(1)./(ytmpl-D_pre)+ABC_l(2)*ytmpl+ABC_l(3))+cc_x_left;
% % % % %     ytmpr=min(floor(D_pre),top_row):-2:btm_row;
% % % % %     xtmpr=-(ABC_r(1)./(ytmpr-D_pre)+ABC_r(2)*ytmpr+ABC_r(3))+cc_x_right;
% % % % % end
% % % % % ytmpl=-(-ytmpl+cc_y+1)/2; %[-1 -240]x640
% % % % % ytmpr=-(-ytmpr+cc_y+1)/2;
% % % % % 
% % % % % ytmpl(round(xtmpl)<=0)=[];
% % % % % xtmpl(round(xtmpl)<=0)=[];
% % % % % ytmpl(round(xtmpl)>ImgSize2)=[];
% % % % % xtmpl(round(xtmpl)>ImgSize2)=[];
% % % % % 
% % % % % ytmpr(round(xtmpr)<=0)=[];
% % % % % xtmpr(round(xtmpr)<=0)=[];
% % % % % ytmpr(round(xtmpr)>ImgSize2)=[];
% % % % % xtmpr(round(xtmpr)>ImgSize2)=[];