%Update particles when no measurement is available
xtmpl=[]; ytmpl=[]; sampled_xtmpl=[]; sampled_ytmpl=[];
xtmpr=[]; ytmpr=[]; sampled_xtmpr=[]; sampled_ytmpr=[];
para=[];
diff_marking_width=0;
logPhai(ii)=Phai0;
logPhaiNaN(ii)=Phai0;
logPhai_f(ii)=Phai0;
logH(ii)=1.665;
logHNaN(ii)=1.665;
logH_f(ii)=1.665;
log_marking_width_left(ii)=log_marking_width_left(ii-1);
log_marking_width_f_left(ii)=log_marking_width_f_left(ii-1);
log_marking_width_right(ii)=log_marking_width_right(ii-1);
log_marking_width_f_right(ii)=log_marking_width_f_right(ii-1);
log_lw(ii)=0;
log_lane_width_f(ii)=lane_width;
lognum_in2(ii)=0;
log_xc_m(ii)=0;
log_d_theta_m(ii)=0;

particle_index=1:P_No;
Particle=Predict_P(:,particle_index);
Estimate_x_mean=mean(Particle(1,:));
Estimate_z_mean=mean(Particle(2,:));
Estimate_theta_mean=mean(Particle(3,:));
[Estimate_xc, Estimate_d_theta]=solve_xc_d_theta(pre_abc1,Estimate_x_mean,Estimate_z_mean,Estimate_theta_mean);
Rtmp=[cos(Estimate_theta_mean-0.5*pi) sin(Estimate_theta_mean-0.5*pi); ...
    -sin(Estimate_theta_mean-0.5*pi) cos(Estimate_theta_mean-0.5*pi)];
%     Particle(1:2,:)=[cos(Estimate_theta_mean-0.5*pi) sin(Estimate_theta_mean-0.5*pi) -Estimate_x_mean; ...
%         -sin(Estimate_theta_mean-0.5*pi) cos(Estimate_theta_mean-0.5*pi) -Estimate_z_mean;]*[Particle(1:2,:); ones(1,P_No)];
Particle(1:2,:)=[Rtmp -Rtmp*[Estimate_x_mean;Estimate_z_mean]]*[Particle(1:2,:); ones(1,P_No)];
Particle(3,:)=Particle(3,:)-(Estimate_theta_mean-0.5*pi);
Particle(4:5,:)=Rtmp*Particle(4:5,:);
log_Estimate_xc(ii)=Estimate_xc;
log_Estimate_d_theta(ii)=Estimate_d_theta;

if leftrightlane==0 && lanelinechange==0
        log_Estimate_xc(ii)=log_Estimate_xc(ii)-lane_width;
%         log_predi_xc(ii)=log_predi_xc(ii)-lane_width;
elseif leftrightlane==2 && lanelinechange==0
%         log_xc_m(ii)=log_xc_m(ii)+lane_width;
        log_Estimate_xc(ii)=log_Estimate_xc(ii)+lane_width;
%         log_predi_xc(ii)=log_predi_xc(ii)+lane_width;
elseif leftrightlane==0 && lanelinechange==1
    leftrightlane=1;
    lanelinechange=0;
elseif leftrightlane==2 && lanelinechange==1
    leftrightlane=1;
    lanelinechange=0;
elseif leftrightlane==1 && lanelinechange==1
    leftrightlane=pre_leftrightlane;
    lanelinechange=0;
    if pre_leftrightlane==0
        log_Estimate_xc(ii)=log_Estimate_xc(ii)-lane_width;
%         log_predi_xc(ii)=log_predi_xc(ii)-lane_width;
    else log_Estimate_xc(ii)=log_Estimate_xc(ii)+lane_width;
%         log_predi_xc(ii)=log_predi_xc(ii)+lane_width;
    end
end
pre_leftrightlane=leftrightlane;

new_pre_World=[Rtmp -Rtmp*[Estimate_x_mean;Estimate_z_mean]]*pre_World([1,3,4],:);
new_pre_World=[new_pre_World(1,:); pre_World(2,:); new_pre_World(2,:); pre_World(4,:)];
new_pre_Xw=new_pre_World(1,:); new_pre_Yw=new_pre_World(2,:); new_pre_Zw=new_pre_World(3,:);
if sum(abs(pre_para(:,2)))==0 %model is hypobola
    Atmp=[new_pre_Zw.^2; new_pre_Zw; ones(1,length(new_pre_Zw))]'; Btmp=new_pre_Xw';
    pre_abc2=pinv(Atmp)*Btmp;
else %model is line
    Atmp=[new_pre_Zw; ones(1,length(new_pre_Zw))]'; Btmp=new_pre_Xw';
    pre_bc=pinv(Atmp)*Btmp;
    pre_abc2=[0; pre_bc];
end
tan_x=Estimate_xc*cos(-Estimate_d_theta);
tan_z=Estimate_xc*sin(-Estimate_d_theta);
Atmp=((new_pre_Zw(end)-tan_z).^2)';
Btmp=(new_pre_Xw(end)-tan_x-tan(Estimate_d_theta)*(new_pre_Zw(end)-tan_z))';
pre_abc1=pinv(Atmp)*Btmp;
pre_abc1=[pre_abc1;tan(Estimate_d_theta)-2*pre_abc1*tan_z; ...
    tan_x-pre_abc1*tan_z^2-(tan(Estimate_d_theta)-2*pre_abc1*tan_z)*tan_z];
World=new_pre_World;
pre_World=World;
continue_frame=continue_frame+1;
% % % 
% % % if leftrightlane==0 && lanelinechange==0
% % %         log_Estimate_xc(ii)=log_Estimate_xc(ii)-lane_width;
% % % %         log_predi_xc(ii)=log_predi_xc(ii)-lane_width;
% % % elseif leftrightlane==2 && lanelinechange==0
% % % %         log_xc_m(ii)=log_xc_m(ii)+lane_width;
% % %         log_Estimate_xc(ii)=log_Estimate_xc(ii)+lane_width;
% % % elseif leftrightlane==0 && lanelinechange==1
% % %     leftrightlane=1;
% % %     lanelinechange=0;
% % % elseif leftrightlane==2 && lanelinechange==1
% % %     leftrightlane=1;
% % %     lanelinechange=0;
% % % elseif leftrightlane==1 && lanelinechange==1
% % %     leftrightlane=pre_leftrightlane;
% % %     lanelinechange=0;
% % %     if pre_leftrightlane==0
% % %         log_Estimate_xc(ii)=log_Estimate_xc(ii)-lane_width;
% % %     else log_Estimate_xc(ii)=log_Estimate_xc(ii)+lane_width;
% % %     end
% % % end
% % % pre_leftrightlane=leftrightlane;

dt=act_t(ii+1);
if v_m~=0
%     v0=v_m+sqrt(var_v).*randn(1,P_No);
    v0=v_m;
else
    v0=0;
end
% phai0=phai_m+sqrt(var_phai).*randn(1,P_No);
phai0=phai_m;
Particle(7,:)=v0;
Particle(8,:)=phai0;
Predict_P(8,:)=Particle(8,:);
Predict_P(7,:)=Particle(7,:);
Predict_P(4,:)=cos(Particle(3,:)).*Particle(7,:);
Predict_P(5,:)=sin(Particle(3,:)).*Particle(7,:);
Predict_P(6,:)=-tan(Particle(8,:)).*Particle(7,:)/L_wheels;
Predict_P(1,:)=Particle(1,:)+Particle(4,:)*dt+0.25*dt/0.2*marking_width*randn(1,P_No)*0;
Predict_P(2,:)=Particle(2,:)+Particle(5,:)*dt+0.25*dt/0.2*marking_width*randn(1,P_No)*0;
Predict_P(3,:)=Particle(3,:)+Particle(6,:)*dt+0.01*dt/0.2*randn(1,P_No)*0;

if sum(Predict_P(2,:)>new_pre_Zw(end))>0.5*length(Predict_P(2,:))
    [Predict_Xc, Predict_d_theta]=C_solve_xc_d_theta(pre_abc2, Predict_P(1:3,:));
else
    [Predict_Xc, Predict_d_theta]=C_solve_xc_d_theta(pre_abc1, Predict_P(1:3,:));
end

predict_x_dot=cos(mean(Particle(3,:)))*v0;
predict_z_dot=sin(mean(Particle(3,:)))*v0;
predict_theta_dot=-tan(phai0)*v0/L_wheels;
predict_x=0+predict_x_dot*dt;
predict_z=0+predict_z_dot*dt;
predict_theta=pi/2+predict_theta_dot*dt;

[est_xtmpl,est_ytmpl,est_xtmpr,est_ytmpr]=Check_next_lane_position(-predict_x,predict_z,pi-predict_theta,Phai,H,pre_abc1,leftrightlane,lane_width, ...
    dist_left_right,L_wheel_cam,Tpix2cam_left, Tpix2cam_right,offset_vector);

logleftrightlane(ii)=leftrightlane;
miss_infor_flag=1;

endt(ii)=toc(start);
endf(ii)=toc(in_start);
