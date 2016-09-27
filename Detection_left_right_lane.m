% determin which lane line the detected line is
sel_xtmpl0=sel_xtmpl;sel_ytmpl0=sel_ytmpl;
xtmpl0=xtmpl;ytmpl0=ytmpl;
sel_xtmpr0=sel_xtmpr;sel_ytmpr0=sel_ytmpr;
xtmpr0=xtmpr;ytmpr0=ytmpr;
num_in0=num_in; para0=para;

if xtmpl0(end)<ImgSize2/2
% % % % %     ir=intersectAB(ytmpl0,sel_ytmpr);
% % % % %     il=intersectAB(sel_ytmpr,ytmpl0);
% % % % %     % toc(start)
% % % % %     sel_ytmpr=sel_ytmpr(ir==1);
% % % % %     sel_xtmpr=sel_xtmpr(ir==1);
% % % % %     Corespx=xtmpl0(il==1);
% % % % %     %         delta_x=round(0.308*(-sel_ytmpr)+7.385);
% % % % %     sel_ytmpr=2*sel_ytmpr-1; %convert back to -480x640 coordinate
% % % % %     ind1=-sel_ytmpr-Offset;
% % % % %     ind1(ind1==0)=1;
% % % % %     delta_x=4*tdx(ind1);
% % % % %     ri_lb=round(sel_xtmpr-delta_x); %right image left boundary
% % % % %     ri_rb=round(sel_xtmpr+delta_x); %right image right boundary
% % % % %     ri_lb(ri_lb<=1)=1;
% % % % %     ri_rb(ri_rb>=(ImgSize2-2))=ImgSize2-2;
% % % % %     li_lb=round(Corespx-2*delta_x);
% % % % %     li_rb=round(Corespx+2*delta_x);
% % % % %     li_lb(li_lb<=1)=1;
% % % % %     li_rb(li_rb>=(ImgSize2-2))=ImgSize2-2;
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
   [~,il,ir]=intersect(ytmpl,sel_ytmpr,'stable');
   sel_ytmpr=sel_ytmpr(ir);
   sel_xtmpr=sel_xtmpr(ir);
   CorepxR=-(sel_xtmpr-cc_x_right);
   CorepxL=-(xtmpl(il)-cc_x_left);
   sel_ytmpr=2*sel_ytmpr-1; %convert back to -480x640 coordinate
    
    disparity=CorepxR-CorepxL;
    indtmp=(disparity~=0);
    Zc=dist_left_right*focal./disparity;
    Xc=CorepxL.*Zc/focal+dist_left_right/2;
    % or Xc=CorepxR.*Zc/focal-dist_left_right/2;
    Yc=(round(sel_ytmpr)+cc_y).*Zc/focal;
    ind_near=(sel_ytmpr<-192); %Use those pixels at 2/3bottom to calculate phai
    
    ind_near=ind_near & indtmp;
    %     if length(Yc)>5
    %         Yc_near=Yc(end-5:end); Zc_near=Zc(end-5:end);
    %     else
    %         Yc_near=Yc; Zc_near=Zc;
    %     end
    if sum(ind_near)>5
        Yc_near=Yc(ind_near); Zc_near=Zc(ind_near); Xc_near=Xc(ind_near);
    else
        Yc_near=Yc(indtmp(end-5:end)); Zc_near=Zc(indtmp(end-5:end)); Xc_near=Xc(indtmp(end-5:end));
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
%     [xtmpl, ytmpl, xtmpr, ytmpr, ~]=refine_model(leftrightlane,para,sel_ytmpr,sel_ytmpl,D_pre,CorepxR,CorepxL);
    
    [~,il,ir]=intersect(ytmpl,ytmpr,'stable');
    CorL=-(xtmpl(il)-cc_x_left);
    CorR=-(xtmpr(ir)-cc_x_right);
    dis=CorR-CorL;
    indtmp=(dis>0);
    
    Zc=dist_left_right*focal./dis;
    Xc=CorL.*Zc/focal+dist_left_right/2;
    Yc=(2*ytmpr(ir)-1+cc_y).*Zc/focal;
    Xc=Xc(indtmp); Yc=Yc(indtmp); Zc=Zc(indtmp);
    
    R=[cos(Theta) -sin(Theta) 0; ...
        cos(Phai)*sin(Theta) cos(Phai)*cos(Theta) -sin(Phai); ...
        sin(Phai)*sin(Theta) sin(Phai)*cos(Theta) cos(Phai);];
    T=[0;H;L_wheel_cam];
    World=[R T+offset_vector;0 0 0 1]*[Xc; Yc; Zc; ones(1,length(Zc))];
    % logXw=World(1); logYw=World(2); logZw=World(3);
    % toc
    Xw=World(1,:); Yw=World(2,:); Zw=World(3,:);
    inf_ind=(isinf(Xw) | isinf(Yw) | isinf(Zw)) | (isnan(Xw) | isnan(Yw) | isnan(Zw)) | (Zw>30);
    Xw(inf_ind)=[];Yw(inf_ind)=[];Zw(inf_ind)=[];
    World(:,inf_ind)=[];
    if sum(paraz(:)==0)
    if sum(abs(para(:,2)))==0 %model is hypobola
        Atmp=[Zw.^2; Zw; ones(1,length(Zw))]';
        Btmp=Xw';
        abc=pinv(Atmp)*Btmp;
        if (abc(1)<0)
            Zwtmp=[Zw(end):-(Zw(end)/length(Zw)):0];
            d_Xw=2*abc(1)*Zw(end)+abc(2);
            const_line=Xw(end)-d_Xw*Zw(end);
            Xwtmp=d_Xw*Zwtmp+const_line;
            Xw=[Xw Xwtmp]; Zw=[Zw Zwtmp];
            Atmp=[Zw.^2; Zw; ones(1,length(Zw))]';
            Btmp=Xw';
            World=[Xw;[Yw zeros(1,length(Xwtmp))];Zw; ones(1,length(Xw))];
            abc=pinv(Atmp)*Btmp;
            Xw=Xw(1:length(Yw));Zw=Zw(1:length(Yw));
            World=World(:,1:length(Yw));
        end
    else %model is line
        Atmp=[Zw; ones(1,length(Zw))]';
        Btmp=Xw';
        bc=pinv(Atmp)*Btmp;
        abc=[0; bc];
    end
    else
        [Xw, Yw, Zw, abc]=equavalent_line(break_location, Xw, Zw);
        World=[Xw;Yw;Zw; ones(1,length(Xw))];
        log_line_para(ii)=3;
    end
    [xc_m, d_theta_m, ~]=solve_xc_d_theta(abc,0,0,pi/2);
else
    % % % % %     il=intersectAB(ytmpr0,sel_ytmpl);
    % % % % %     ir=intersectAB(sel_ytmpl,ytmpr0);
    % % % % %     % toc(start)
    % % % % %     sel_ytmpl=sel_ytmpl(il==1);
    % % % % %     sel_xtmpl=sel_xtmpl(il==1);
    % % % % %     Corespx=xtmpr0(ir==1);
    % % % % %     %         delta_x=round(0.308*(-sel_ytmpl)+7.385);
    % % % % %     sel_ytmpl=2*sel_ytmpl-1; %convert back to -480x640 coordinate
    % % % % %     ind1=-sel_ytmpl-Offset;
    % % % % %     ind1(ind1==0)=1;
    % % % % %     delta_x=4*tdx(ind1);
    % % % % %     li_lb=round(sel_xtmpl-delta_x); %right image left boundary
    % % % % %     li_rb=round(sel_xtmpl+delta_x); %right image right boundary
    % % % % %     li_lb(li_lb<=1)=1;
    % % % % %     li_rb(li_rb>=(ImgSize2-2))=ImgSize2-2;
    % % % % %     ri_lb=round(Corespx-2*delta_x);
    % % % % %     ri_rb=round(Corespx+2*delta_x);
    % % % % %     ri_lb(ri_lb<=1)=1;
    % % % % %     ri_rb(ri_rb>=(ImgSize2-2))=ImgSize2-2;
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
    [~,ir,il]=intersect(ytmpr,sel_ytmpl,'stable');
    sel_ytmpl=sel_ytmpl(il);
    sel_xtmpl=sel_xtmpl(il);
    CorepxL=-(sel_xtmpl-cc_x_left);
    CorepxR=-(xtmpr(ir)-cc_x_right);
    sel_ytmpl=2*sel_ytmpl-1; %convert back to -480x640 coordinate
    
    disparity=CorepxR-CorepxL;
    indtmp=(disparity~=0);
    Zc=dist_left_right*focal./disparity;
    Xc=CorepxL.*Zc/focal+dist_left_right/2;
    % or Xc=CorepxR.*Zc/focal-dist_left_right/2;
    Yc=(round(sel_ytmpl)+cc_y).*Zc/focal;
    ind_near=(sel_ytmpl<-192); %Use those pixels at 2/3bottom to calculate phai
    ind_near=ind_near & indtmp;
    %     if length(Yc)>5
    %         Yc_near=Yc(end-5:end); Zc_near=Zc(end-5:end);
    %     else
    %         Yc_near=Yc; Zc_near=Zc;
    %     end
    if sum(ind_near)>5
        Yc_near=Yc(ind_near); Zc_near=Zc(ind_near); Xc_near=Xc(ind_near);
    else
        Yc_near=Yc(indtmp(end-5:end)); Zc_near=Zc(indtmp(end-5:end)); 
        Xc_near=Xc(indtmp(end-5:end));
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
%     [xtmpl, ytmpl, xtmpr, ytmpr, ~]=refine_model(leftrightlane,para,sel_ytmpr,sel_ytmpl,D_pre,CorepxR,CorepxL);
    
    [~,il,ir]=intersect(ytmpl,ytmpr,'stable');
    CorL=-(xtmpl(il)-cc_x_left);
    CorR=-(xtmpr(ir)-cc_x_right);
    dis=CorR-CorL;
    indtmp=(dis>0);
    
    Zc=dist_left_right*focal./dis;
    Xc=CorL.*Zc/focal+dist_left_right/2;
    Yc=(2*ytmpl(il)-1+cc_y).*Zc/focal;
    Xc=Xc(indtmp); Yc=Yc(indtmp); Zc=Zc(indtmp);
    
    R=[cos(Theta) -sin(Theta) 0; ...
        cos(Phai)*sin(Theta) cos(Phai)*cos(Theta) -sin(Phai); ...
        sin(Phai)*sin(Theta) sin(Phai)*cos(Theta) cos(Phai);];
    T=[0;H;L_wheel_cam];
    World=[R T+offset_vector;0 0 0 1]*[Xc; Yc; Zc; ones(1,length(Zc))];
    % logXw=World(1); logYw=World(2); logZw=World(3);
    % toc
    Xw=World(1,:); Yw=World(2,:); Zw=World(3,:);
    inf_ind=(isinf(Xw) | isinf(Yw) | isinf(Zw)) | (isnan(Xw) | isnan(Yw) | isnan(Zw)) | (Zw>30);
    Xw(inf_ind)=[];Yw(inf_ind)=[];Zw(inf_ind)=[];
    World(:,inf_ind)=[];
    if sum(paraz(:)==0)
    if sum(abs(para(:,2)))==0 %model is hypobola
        Atmp=[Zw.^2; Zw; ones(1,length(Zw))]';
        Btmp=Xw';
        abc=pinv(Atmp)*Btmp;
        if (abc(1)>0)
            Zwtmp=[Zw(end):-(Zw(end)/length(Zw)):0];
            d_Xw=2*abc(1)*Zw(end)+abc(2);
            const_line=Xw(end)-d_Xw*Zw(end);
            Xwtmp=d_Xw*Zwtmp+const_line;
            Xw=[Xw Xwtmp]; Zw=[Zw Zwtmp];
            Atmp=[Zw.^2; Zw; ones(1,length(Zw))]';
            Btmp=Xw';
            World=[Xw;[Yw zeros(1,length(Xwtmp))];Zw; ones(1,length(Xw))];
            abc=pinv(Atmp)*Btmp;
            Xw=Xw(1:length(Yw));Zw=Zw(1:length(Yw));
            World=World(:,1:length(Yw));
        end
    else %model is line
        Atmp=[Zw; ones(1,length(Zw))]';
        Btmp=Xw';
        bc=pinv(Atmp)*Btmp;
        abc=[0; bc];
    end
    else
        [Xw, Yw, Zw, abc]=equavalent_line(break_location, Xw, Zw);
       World=[Xw;Yw;Zw; ones(1,length(Xw))];
       log_line_para(ii)=3;
    end
    [xc_m, d_theta_m, ~]=solve_xc_d_theta(abc,0,0,pi/2);
end
% % % % consider 3 lane lines
% % % xc_m1=xc_m-lane_width; xc_m2=xc_m; xc_m3=xc_m+lane_width;
% % % del_dist_last_xc=abs([xc_m1 xc_m2 xc_m3]-log_Estimate_xc(ii-1));
% % % [~,i_min]=min(del_dist_last_xc);
% % % if i_min==1
% % %     leftrightlane=0;
% % % elseif i_min==2
% % %     leftrightlane=1;
% % % else leftrightlane=2;
% % % end

% % % % % consider 2 lane lines only
xc_m1=xc_m-lane_width; xc_m2=xc_m;
del_dist_last_xc=abs([xc_m1 xc_m2]-log_Estimate_xc(ii-1));
[~,i_min]=min(del_dist_last_xc);
if i_min==1
    leftrightlane=0;
elseif i_min==2
    leftrightlane=1;
end

sel_xtmpl=sel_xtmpl0;sel_ytmpl=sel_ytmpl0;
xtmpl=xtmpl0;ytmpl=ytmpl0;
sel_xtmpr=sel_xtmpr0;sel_ytmpr=sel_ytmpr0;
xtmpr=xtmpr0;ytmpr=ytmpr0;
num_in=num_in0; para=para0;