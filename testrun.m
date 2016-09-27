clear all;
warning off;

% leftimg_dir='D:\Stereo\Image Collection4\left';
% rightimg_dir='D:\Stereo\Image Collection4\right';
leftimg_dir='D:\Stereo\Image Collection6\right';
rightimg_dir='D:\Stereo\Image Collection6\left';
save_dir='D:\Stereo\Image Collection6';

global_varibles;
InitializationGlobalVariable;
load 'cali3.mat';
%matlabpool(4);
D_pre=852*tan(15/180*pi);H=1.665;marking_width_pre=0.1;diff_marking_width=0;
marking_width_pre2=0.1;L_wheel_cam=0.675; lane_width=3.0;
top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even
cc_y=2*ceil(cc_y/2);
% D_pre=174.4855; H=1.6428;
SampleRectification;
Ridge_width_ini=12.5;
ImgSize1=480; ImgSize2=640;
RidgeParameter;
wd=zeros(5,5,8); wd(3,3,:)=1;
wd(3,:,1)=1;
wd(4,1:2,2)=1; wd(2,4:5,2)=1;wd(3,2,2)=1;wd(3,4,2)=1;
wd(1,5,3)=1;wd(2,4,3)=1;wd(4,2,3)=1;wd(5,1,3)=1;
wd(4:5,2,4)=1;wd(1:2,4,4)=1;wd(4,3,4)=1;wd(2,3,4)=1;
wd(:,3,5)=1;
wd(4:5,4,6)=1;wd(1:2,2,6)=1;wd(4,3,6)=1;wd(2,3,6)=1;
wd(1,1,7)=1;wd(2,2,7)=1;wd(4,4,7)=1;wd(5,5,7)=1;
wd(2,1:2,8)=1; wd(4,4:5,8)=1; wd(3,2,8)=1; wd(3,4,8)=1;
convimg=zeros(240,1310,8);
kcsditmpl=zeros(ImgSize1/2,ImgSize2);
kcsditmpr=zeros(ImgSize1/2,ImgSize2);
xtmpl=[]; xtmpr=[];
% ite_s=1712; ite_e=1800;
% ite_s=52; ite_e=405;
% ite_s=646; ite_e=845;
% ite_s=961; ite_e=1255;
% ite_s=1516; ite_e=1800;
% Theta=-5/180*pi; %for Image Collection 1 and 2

% Image Collection 3
% ite_s=516; ite_e=560;
% ite_s=571; ite_e=930;
% ite_s=1041; ite_e=1230;
% ite_s=1321; ite_e=1560;
% ite_s=1681; ite_e=2050;
% ite_s=2166; ite_e=2460;

%Image Collection 4
% ite_s=1465; ite_e=2153;
% ite_s=1; ite_e=400;
% ite_s=516; ite_e=965;
ite_s=1080; ite_e=1345;
% ite_s=1435; ite_e=1700;
% ite_s=1834; ite_e=2250;
% ite_s=2368; ite_e=2680;
Theta=-0.8/180*pi; %for Image Collection 3,4,5

for ii=ite_s:ite_e
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
    
    ImgSmooth=zeros(size(img));
    for i=1:(length(ix)-1)
        ImgSmooth(ix(i):(ix(i+1)-1),:)=conv2(img(ix(i):(ix(i+1)-1),:),xfil2{i},'same');
    end
    
% % %     if diff_marking_width==1
% % %         for i=start_ix2:(length(ix2)-1)
% % %             ImgSmooth(ix2(i):(ix2(i+1)-1),[x_conv_li2_lb(i):x_conv_li2_rb(i),(x_conv_ri2_lb(i)+ImgSize2+30):(x_conv_ri2_rb(i)+ImgSize2+30)])= ...
% % %                 conv2(img(ix2(i):(ix2(i+1)-1),[x_conv_li2_lb(i):x_conv_li2_rb(i),(x_conv_ri2_lb(i)+ImgSize2+30):(x_conv_ri2_rb(i)+ImgSize2+30)]),xfil2_2{i},'same');
% % %         end
% % %     end
    
    ImgSmooth=conv2(ImgSmooth,yfil2,'same');
    [kcsditmp,th]=C_Ridge(ImgSmooth,306400);
    for i=1:8
        convimg(:,:,i)=conv2(kcsditmp,wd(:,:,i),'same');
    end
    maxconvimg=max(convimg,[],3);

    dumy1=maxconvimg;
    dumy1(dumy1<4.5*th)=0; dumy1(dumy1>=4.5*th)=1;
    % tic;
    num_ite=15000; %number of iterations for RANSAC
    if ~isempty(xtmpl) && ~isempty(xtmpr)
        maxconvimg=C_Search_region(maxconvimg,xtmpl,xtmpr+ImgSize2+30,-ytmpr(1),-ytmpr(end));
        maxconvimg(-ytmpr(end):end,:)=0;
        maxconvimg(1:-ytmpr(1),:)=0;
        num_ite=960;
    end
    % toc
    dumy2=maxconvimg;
    dumy2(dumy2<4.5*th)=0;dumy2(dumy2>=4.5*th)=1;
    dumy=[dumy1;dumy2];
    % dumy(:,ImgSize2/2:ImgSize2)=0;
    % dumy(:,end-320:end)=0;
    imwrite(dumy,strcat(save_dir,'\Ridge',int2str(ii),'.jpg'),'jpg');
    
    [sel_xtmpl,sel_ytmpl, xtmpl, ytmpl,sel_xtmpr,sel_ytmpr, xtmpr, ytmpr, num_in, para]=linefitting_simu(maxconvimg(:,1:ImgSize2),maxconvimg(:,end-ImgSize2+1:end),th, D_pre, H, num_ite);
    lognum_in(ii)=num_in;
    if num_in<100
        diff_marking_width=0;
        logPhai(ii)=0.2757;
        logH(ii)=1.665;
        log_marking_width(ii)=marking_width_pre;
        log_marking_width_f(ii)=marking_width_pre;
        log_lane_width_f(ii)=0;
        endt(ii)=toc(start);
        endf(ii)=toc(in_start);
        continue;
    end
    %convert back to original image coordinate 480x640
    % % [~,ir,il]=intersect(sel_ytmpr, ytmpl, 'stable');
    if mean(xtmpr)<320 %This condition may not be true...Need to update this based on Xworld coordinate
        leftrightlane=0;
    else
        leftrightlane=1;
    end
    if leftrightlane==0 %left lane
        ir=intersectAB(ytmpl,sel_ytmpr);
        il=intersectAB(sel_ytmpr,ytmpl);
        % toc(start)
        sel_ytmpr=sel_ytmpr(ir==1);
        sel_xtmpr=sel_xtmpr(ir==1);
        Corespx=xtmpl(il==1);
        delta_x=round(0.308*(-sel_ytmpr)+7.385);
        sel_ytmpr=2*sel_ytmpr-1; %convert back to -480x640 coordinate
        ri_lb=round(sel_xtmpr-delta_x); %right image left boundary
        ri_rb=round(sel_xtmpr+delta_x); %right image right boundary
        ri_lb(ri_lb<=1)=1;
        ri_rb(ri_rb>=ImgSize2)=ImgSize2;
        li_lb=round(Corespx-2*delta_x);
        li_rb=round(Corespx+2*delta_x);
        li_lb(li_lb<=1)=1;
        li_rb(li_rb>=ImgSize2)=ImgSize2;
        ri_lb=ri_lb-1;
        ri_rb=ri_rb-1;
        li_lb=li_lb-1;
        li_rb=li_rb-1;
        
        IntensityR=rightimg(-sel_ytmpr,:);
        CorespL=leftimg(-sel_ytmpr,:);
        % toc(start)
        disparitytmp=ZNCC(ri_lb', ri_rb', li_lb', li_rb', IntensityR, CorespL);
        % toc(start)
        li_lb=li_lb+1;
        ri_lb=ri_lb+1;
        % Refer to the sign convention in 3D Vision notes
        % Convert to left and right camera coordinates respectively
        CorepxL=-(li_lb+disparitytmp'+(round(sel_xtmpr)-ri_lb)-cc_x_left);
        CorepxR=-(round(sel_xtmpr)-cc_x_right);
    else %right lane
        il=intersectAB(ytmpr,sel_ytmpl);
        ir=intersectAB(sel_ytmpl,ytmpr);
        % toc(start)
        sel_ytmpl=sel_ytmpl(il==1);
        sel_xtmpl=sel_xtmpl(il==1);
        Corespx=xtmpr(ir==1);
        delta_x=round(0.308*(-sel_ytmpl)+7.385);
        sel_ytmpl=2*sel_ytmpl-1; %convert back to -480x640 coordinate
        li_lb=round(sel_xtmpl-delta_x); %right image left boundary
        li_rb=round(sel_xtmpl+delta_x); %right image right boundary
        li_lb(li_lb<=1)=1;
        li_rb(li_rb>=ImgSize2)=ImgSize2;
        ri_lb=round(Corespx-2*delta_x);
        ri_rb=round(Corespx+2*delta_x);
        ri_lb(ri_lb<=1)=1;
        ri_rb(ri_rb>=ImgSize2)=ImgSize2;
        li_lb=li_lb-1;
        li_rb=li_rb-1;
        ri_lb=ri_lb-1;
        ri_rb=ri_rb-1;
        
        IntensityL=leftimg(-sel_ytmpl,:);
        CorespR=rightimg(-sel_ytmpl,:);
        % toc(start)
        disparitytmp=ZNCC(li_lb', li_rb', ri_lb', ri_rb', IntensityL, CorespR);
        % toc(start)
        li_lb=li_lb+1;
        ri_lb=ri_lb+1;
        % Refer to the sign convention in 3D Vision notes
        % Convert to left and right camera coordinates respectively
        CorepxR=-(ri_lb+disparitytmp'+(round(sel_xtmpl)-li_lb)-cc_x_right);
        CorepxL=-(round(sel_xtmpl)-cc_x_left);
    end
    disparity=CorepxR-CorepxL;
    indtmp=(disparity~=0);
    Zc=dist_left_right*focal./disparity;
    Xc=CorepxL.*Zc/focal+dist_left_right/2;
    % or Xc=CorepxR.*Zc/focal-dist_left_right/2;
    if leftrightlane==0
        Yc=(round(sel_ytmpr)+cc_y).*Zc/focal;
        ind_near=(sel_ytmpr<-192); %Use those pixels at 3/5bottom to calculate phai
    else
        Yc=(round(sel_ytmpl)+cc_y).*Zc/focal;
        ind_near=(sel_ytmpl<-192); %Use those pixels at 3/5bottom to calculate phai
    end
    ind_near=ind_near & indtmp;
    %     if length(Yc)>5
    %         Yc_near=Yc(end-5:end); Zc_near=Zc(end-5:end);
    %     else
    %         Yc_near=Yc; Zc_near=Zc;
    %     end
    if sum(ind_near)>5
        Yc_near=Yc(ind_near); Zc_near=Zc(ind_near); Xc_near=Xc(ind_near);
    else
        Yc_near=Yc(indtmp); Zc_near=Zc(indtmp); Xc_near=Xc(indtmp);
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
    [xtmpl, ytmpl, xtmpr, ytmpr, refine_para]=refine_model(leftrightlane,para,sel_ytmpr,sel_ytmpl,D_pre,CorepxR,CorepxL);
    if numel(ytmpl)==0 || numel(ytmpr)==0
        endt(ii)=toc(start);
        endf(ii)=toc(in_start);
        continue;
    end
    [sel_xtmpl2,sel_ytmpl2,xtmpl2, ytmpl2, sel_xtmpr2,sel_ytmpr2,xtmpr2, ytmpr2, num_in2, para2, lwtmp]=fit_other_lane(leftrightlane,dumy1,refine_para,-ytmpr(1),-ytmpr(end),H,Phai,D_pre);
    lognum_in2(ii)=num_in2;
    if isempty(lwtmp)
        log_lw(ii,:)=0;
    else
        log_lw(ii,:)=lwtmp;
    end
    log_non_0_lw=log_lw(log_lw~=0);
    log_non_0_num_in2=lognum_in2(log_lw~=0);
    m=-min([10, length(log_non_0_lw)]):-1;
    coe_fil=cos(m*pi/2/max(-m));
    coe_fil=log_non_0_num_in2((end-max(-m)+1):end).*coe_fil;
    coe_fil=coe_fil'/sum(coe_fil);
    lane_width=sum(coe_fil.*log_non_0_lw((end-max(-m)+1):end));
    log_lane_width_f(ii)=lane_width;
    
    R=[cos(Theta) -sin(Theta) 0; ...
        cos(Phai)*sin(Theta) cos(Phai)*cos(Theta) -sin(Phai); ...
        sin(Phai)*sin(Theta) sin(Phai)*cos(Theta) cos(Phai);];
    T=[0;H;L_wheel_cam]; %0.675m is the horizontal distance from the vehicle control point to the camera
    World=[R T;0 0 0 1]*[Xc; Yc; Zc; ones(1,length(Zc))];
    Xw=World(1,:); Yw=World(2,:); Zw=World(3,:);
    inf_ind=(isinf(Xw) | isinf(Yw) | isinf(Zw)) | (isnan(Xw) | isnan(Yw) | isnan(Zw));
    Xw(inf_ind)=[];Yw(inf_ind)=[];Zw(inf_ind)=[];
    Atmp=[Zw.^2; Zw; ones(1,length(Zw))]'; Btmp=Xw';
    abc=pinv(Atmp)*Btmp;
    [xc, d_theta]=solve_xc_d_theta(abc,0,0,pi/2);
    log_xc(ii)=xc;
    log_d_theta(ii)=d_theta;
    

    % logXw=World(1); logYw=World(2); logZw=World(3);
    % toc
    
    if num_in2>num_in
        xtmpl=xtmpl2; ytmpl=ytmpl2; xtmpr=xtmpr2; ytmpr=ytmpr2;
    end
    %     logWidth(ii)=f_width;
    
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
        log_marking_width(ii)=marking_width;
        marking_width=median(log_marking_width(max([ite_s, ii-4]):ii));
        log_marking_width_f(ii)=marking_width;
        if marking_width~=marking_width_pre
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
        end
    else
        log_marking_width(ii)=marking_width_pre;
        log_marking_width_f(ii)=marking_width_pre;
    end
    
% % %     if num_in2>50
% % %         if leftrightlane==0
% % %             ind_near2=(sel_ytmpl2<-160); %Use those pixels at 2/3bottom to calculate phai
% % %         else
% % %             ind_near2=(sel_ytmpr2<-160); %Use those pixels at 2/3bottom to calculate phai
% % %         end
% % %         if sum(ind_near2)>6
% % %             ind_near2(1:(end-7))=0;
% % %         end
% % %         if mod(sum(ind_near2),2)==0
% % %             indtmp=find(ind_near2==1,1);
% % %             ind_near2(indtmp)=0;
% % %         end
% % %         if leftrightlane==0
% % %             y_for_tdx2=sel_ytmpl2(ind_near2)';
% % %             x_for_tdx2=round(sel_xtmpl2(ind_near2)');
% % %             Intensity2=rightimg(-y_for_tdx2,:);
% % %         else
% % %             y_for_tdx2=sel_ytmpr2(ind_near2)';
% % %             x_for_tdx2=round(sel_xtmpr2(ind_near2)');
% % %             Intensity2=leftimg(-y_for_tdx2,:);
% % %         end
% % %         y_for_tdx2=y_for_tdx2+256; %convert to camera coordinate
% % %         tdxtmp1=round((cos(Phai)/H*(D_pre-y_for_tdx2)*0.10)/2);
% % %         tdxtmp2=round((cos(Phai)/H*(D_pre-y_for_tdx2)*0.15)/2);
% % %         tdxtmp3=round((cos(Phai)/H*(D_pre-y_for_tdx2)*0.20)/2);
% % %         quater_tdxtmp3=round(1/4*tdxtmp3);
% % %         x_for_tdx_lb=x_for_tdx2-tdxtmp3-quater_tdxtmp3;
% % %         x_for_tdx_rb=x_for_tdx2+tdxtmp3+quater_tdxtmp3;
% % %         x_for_tdx_lb(x_for_tdx_lb<1)=1;
% % %         x_for_tdx_lb(x_for_tdx_lb>ImgSize2)=ImgSize2;
% % %         x_for_tdx_rb(x_for_tdx_rb<1)=1;
% % %         x_for_tdx_rb(x_for_tdx_rb>ImgSize2)=ImgSize2;
% % %         marking_width_all2=C_find_lane_marking_width(Intensity2,x_for_tdx_lb,x_for_tdx_rb,tdxtmp1,tdxtmp2,tdxtmp3,x_for_tdx2);
% % %         marking_width2=median(marking_width_all2);
% % %         log_marking_width2(ii)=marking_width2;
% % %         marking_width2=median(log_marking_width2(max([ite_s, ii-4]):ii));
% % %         log_marking_width_f2(ii)=marking_width2;
% % %         if marking_width2~=marking_width_pre2;
% % %             marking_width_pre2=marking_width2;
% % %             if mod(floor(D_pre),2)==0
% % %                 D_pre1=floor(D_pre)-1;
% % %             else D_pre1=floor(D_pre);
% % %             end
% % %             yy=min(floor(D_pre1),215):-1:-224;
% % %             tdxtmp=(cos(Phai)/H*(D_pre-yy)*marking_width2)/2; %0.10m is the lane marking width
% % %             if length(tdxtmp)>=length(tdx)
% % %                 tdx=round(tdxtmp((end-length(tdx)+1):end));
% % %             else
% % %                 tdx=round([0.5*ones(1,(length(tdx)-length(tdxtmp))) tdxtmp]);
% % %             end
% % %             [~,ix2,~]=unique(tdx,'stable');ix2=ix2';
% % %             j=1;
% % %             for i=ix2
% % %                 weight3=exp((-0.5*(-tdx(i):tdx(i)).^2/(tdx(i))^2));
% % %                 weight3=weight3/sum(weight3);
% % %                 xfil2_2{j}=weight3;
% % %                 j=j+1;
% % %             end
% % %             %     ix2((-ix2+256)>D_pre1)=[];
% % %             ix2(length(ix2)+1)=ImgSize1-Offset+1;
% % %             ix2=round(ix2/2)+Offset/2;
% % %             start_ix2=find((-ix2*2+256)<D_pre1,1);
% % %             if para2(:,2)==0 %model is hypobola
% % %                 x_conv_li2=para2(1)./((-ix2*2+256)-D_pre)+(-ix2*2+256)*para2(2)+para2(3);
% % %                 x_conv_ri2=para2(1)./((-ix2*2+256)-D_pre)+(-ix2*2+256)*para2(4)+para2(5);
% % %             else %model is line
% % %                 x_conv_li2=-para2(2,2)*(-ix2*2+256)-0.5*para2(3,2);
% % %                 x_conv_ri2=-para2(4,2)*(-ix2*2+256)-0.5*para2(5,2);
% % %             end
% % %             x_conv_li2=-x_conv_li2+cc_x_left; x_conv_ri2=-x_conv_ri2+cc_x_right;
% % %             x_conv_li2_lbtmp=round(x_conv_li2-81); x_conv_li2_lbtmp(x_conv_li2_lbtmp<1)=1;x_conv_li2_lbtmp(x_conv_li2_lbtmp>ImgSize2)=ImgSize2;
% % %             x_conv_li2_rbtmp=round(x_conv_li2+81); x_conv_li2_rbtmp(x_conv_li2_rbtmp<1)=1;x_conv_li2_rbtmp(x_conv_li2_rbtmp>ImgSize2)=ImgSize2;
% % %             x_conv_ri2_lbtmp=round(x_conv_ri2-81); x_conv_ri2_lbtmp(x_conv_ri2_lbtmp<1)=1;x_conv_ri2_lbtmp(x_conv_ri2_lbtmp>ImgSize2)=ImgSize2;
% % %             x_conv_ri2_rbtmp=round(x_conv_ri2+81); x_conv_ri2_rbtmp(x_conv_ri2_rbtmp<1)=1;x_conv_ri2_rbtmp(x_conv_ri2_rbtmp>ImgSize2)=ImgSize2;
% % %             
% % %             delta=x_conv_li2_lbtmp(1:(end-1))-x_conv_li2_lbtmp(2:end);
% % %             indt=find(delta<0);
% % %             x_conv_li2_lb=zeros(length(delta),1);
% % %             x_conv_li2_lb(indt)=x_conv_li2_lbtmp(indt);
% % %             indt=find(delta>=0);
% % %             x_conv_li2_lb(indt)=x_conv_li2_lbtmp(indt+1);
% % %             
% % %             delta=x_conv_li2_rbtmp(1:(end-1))-x_conv_li2_rbtmp(2:end);
% % %             indt=find(delta>=0);
% % %             x_conv_li2_rb=zeros(length(delta),1);
% % %             x_conv_li2_rb(indt)=x_conv_li2_rbtmp(indt);
% % %             indt=find(delta<0);
% % %             x_conv_li2_rb(indt)=x_conv_li2_rbtmp(indt+1);
% % %             
% % %             delta=x_conv_ri2_lbtmp(1:(end-1))-x_conv_ri2_lbtmp(2:end);
% % %             indt=find(delta<0);
% % %             x_conv_ri2_lb=zeros(length(delta),1);
% % %             x_conv_ri2_lb(indt)=x_conv_ri2_lbtmp(indt);
% % %             indt=find(delta>=0);
% % %             x_conv_ri2_lb(indt)=x_conv_ri2_lbtmp(indt+1);
% % %             
% % %             delta=x_conv_ri2_rbtmp(1:(end-1))-x_conv_ri2_rbtmp(2:end);
% % %             indt=find(delta>=0);
% % %             x_conv_ri2_rb=zeros(length(delta),1);
% % %             x_conv_ri2_rb(indt)=x_conv_ri2_rbtmp(indt);
% % %             indt=find(delta<0);
% % %             x_conv_ri2_rb(indt)=x_conv_ri2_rbtmp(indt+1);
% % %         end
% % %     else
% % %         log_marking_width2(ii)=marking_width_pre2;
% % %         marking_width_pre2=marking_width_pre;
% % %     end
 
% % %     % for the other lane converlution:
% % %     if marking_width_pre~=marking_width_pre2
% % %         diff_marking_width=1;
% % %     else
% % %         diff_marking_width=0;
% % %     end
    endt(ii)=toc(start);
    endf(ii)=toc(in_start);
    plotimg;
end
