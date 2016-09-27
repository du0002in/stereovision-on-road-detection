clear all;pause(0.01);
warning off;
f='D:\StereoLaptop_new\Image Collection(new Cam)8';
save_dir=f;
% f='D:\StereoLaptop_new\RealRun_Laptop_new_Cam\RealRunImg\test166';
leftimg_dir=[f '\left'];
rightimg_dir=[f '\right'];
% save_dir=[f '\Offline'];
% % % leftimg_dir='D:\Stereo\RealRun_Laptop_new_Cam\RealRunImg\test93\left';
% % % rightimg_dir='D:\Stereo\RealRun_Laptop_new_Cam\RealRunImg\test93\right';
% % % save_dir='D:\Stereo\RealRun_Laptop_new_Cam\RealRunImg\test93\Offline';
% leftimg_dir='D:\Stereo\RealRun_Laptop\RealRunImg\test15\left';
% rightimg_dir='D:\Stereo\RealRun_Laptop\RealRunImg\test15\right';
% save_dir='D:\Stereo\RealRun_Laptop\RealRunImg\test15\offline';
% leftimg_dir='D:\Stereo\MPC\GA with Simple Map\Simulated Img\left';
% rightimg_dir='D:\Stereo\MPC\GA with Simple Map\Simulated Img\right';
% save_dir='D:\Stereo\MPC\GA with Simple Map\Simulated Img';
% leftimg_dir='D:\Stereo\Calibration\realcali\realleft';
% rightimg_dir='D:\Stereo\Calibration\realcali\realright';
% save_dir='D:\Stereo\Calibration\realcali';
figure;
%for Image Collection(new Cam)1
log_ite_s=[3 193 461 3];
log_ite_e=[378 857 857 911];

%for Image Collection(new Cam)2,3,4
% log_ite_s=[1    170 318 459 622 770];
% log_ite_e=[113  250 381 541 710 860];
load 'pca_temp3.mat';
load 'pca_temp3_arrow_tip2.mat';
% % load 'pca_temp3_line.mat';
load 'True_ind.mat';
load 'ROI2.mat';
for outloop=1:1
    ite_s=log_ite_s(outloop);
    ite_e=log_ite_e(outloop);
    global_varibles;
    InitializationGlobalVariable;
    load 'cali3.mat'; %'cali3_old2.mat' is for img collection 6. The original webcam calibration results
    load(strcat(save_dir,'\odemotry.mat'));%actural data readings from EV velocity and steering
    % load 'fft_boundary.mat';
    save_dir=f;
    miss_infor_flag=0;
    Phai0=0.3070;
    D_pre=focal*tan(Phai0);H=1.665;marking_width_pre=0.10;diff_marking_width=0;lane_width=2.95;
    marking_width_pre2=0.10;continue_frame=0;marking_width=0.1;
    top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
    btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even
    cc_y=2*ceil(cc_y/2);
    % D_pre=174.4855; H=1.6428;
    SampleRectification;
    Ridge_width_ini=8.5;
    ImgSize1=480; ImgSize2=640; NFFT=2^nextpow2(ImgSize2);
    fre_seq=ImgSize2/2*linspace(0,1,NFFT/2+1);
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
    %Image Collection 2
    % ite_s=52; ite_e=405;
    % ite_s=646; ite_e=845;
    % ite_s=961; ite_e=1255;
    % ite_s=1516; ite_e=1800;
    % Theta=-5/180*pi;
    % This block is for template matching initialization
    load 'roadmarkingtemplates_low.mat';
    del_w=-30:5:30;
    theta=asin(del_w/65);
    del_h=floor(26*cos(theta));
    for i=1:length(del_w)
        log_ind{i}=[round((65*sin(theta(i))-1)/(26*cos(theta(i))-1)*((1:del_h(i))-1)+1); 1:del_h(i)];
        x=[log_ind{i}(1,:) log_ind{i}(1,:)+1 log_ind{i}(1,:)+2];
        y=[log_ind{i}(2,:) log_ind{i}(2,:) log_ind{i}(2,:)];
        min_x=min(x);
        if min_x<=0
            x=x-min_x+1;
        end
        mask=zeros(max(y), max(x));
        ind=sub2ind(size(mask),y,x);
        mask(ind)=1;
        log_mask{i}=mask;
    end
    log_atrribute=cell(1,30);
    [ximg,yimg]=meshgrid(1:640,1:2:480);
    xxc=-ximg+cc_x_right;
    yyc=-yimg+cc_y;
    total_counter=0;
    
    % Image Collection 3
    % ite_s=1; ite_e=415;
    % ite_s=571; ite_e=930;
    % ite_s=1041; ite_e=1230;
    % ite_s=1321; ite_e=1560;
    % ite_s=1681; ite_e=2050;
    % ite_s=2166; ite_e=2460;
    
    %Image Collection 4
    % ite_s=516; ite_e=560;
    % ite_s=2150; ite_e=2220;
    % ite_s=1; ite_e=400;
    % ite_s=516; ite_e=965;
    % ite_s=1100; ite_e=1345;
    % ite_s=1435; ite_e=1700;
    % ite_s=1834; ite_e=2250;
    % ite_s=2368; ite_e=2680;
    Theta=-0.002;
    % offset_vector=[-0.09;0.01;-0.34]; %World=R*Cam+offset_vector;
    offset_vector=[-0.0;0.01;-0.34]; %World=R*Cam+offset_vector;
    [Tpix2cam_left, Tpix2cam_right]=Pix2Cam(focal,cc_x_left,cc_x_right,cc_y);
    
    pre_line_type=1;
    % ite_s=1384; ite_e=1498;
    % vid1 = videoinput('winvideo', 1, 'YUY2_640x480');
    % set(vid1,'FramesPerTrigger',1);
    % set(vid1,'TriggerRepeat',Inf);
    % triggerconfig(vid1,'manual');
    % start(vid1);
    % preview(vid1);
    % log_dumy=zeros(1,2);
    ximg=[]; yimg=[]; log_xyimg{1}=[];
    %     try
    KickOff_PF2;
    % intesnity_td_min=230;
    % intesnity_td_max=255;
    log_p_tip(ii)=0;
    for ii=(ite_s+1):ite_e
        % for ii=(93):ite_e
        if ii==8
            a=0;
        end
        in_start=tic;
        %     trigger(vid1);
        %     dd=getdata(vid1);
        % [v_m, phai_m]=get_measurements;
        v_m=act_v(ii); phai_m=act_phai(ii);
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
        mask=img;mask(img<intesnity_td_min)=0;mask(img>=intesnity_td_min)=1;mask(img>intesnity_td_max)=0;
        %     C_in_region(dumy,xtmpli_lb,xtmpli_rb,xtmpri_lb+ImgSize2+30,xtmpri_rb+ImgSize2+30,ytopr,ybtmr);
        if ~isempty(log_expos_box{ii-1})
            expos_box=log_expos_box{ii-1};
        end
        expo_mask=C_in_region(img, expos_box{1}+40, expos_box{2}-40, expos_box{3}+ImgSize2+30+40, expos_box{4}+ImgSize2+30-40, expos_box{5}, expos_box{6});
        expo_mask(1:25,:)=0;
        expo_mask_l1=expo_mask(25:55,1:ImgSize2);
        expo_mask_l2=expo_mask(56:120,1:ImgSize2);
        expo_mask_l3=expo_mask(121:end,1:ImgSize2);
        log_road_bright_l(ii)=mean(expo_mask_l1(expo_mask_l1~=0))+ ...
            mean(expo_mask_l2(expo_mask_l2~=0))+ ...
            mean(expo_mask_l3(expo_mask_l3~=0));
        log_road_bright_l(ii)=log_road_bright_l(ii)/3;
        
        expo_mask_r1=expo_mask(25:55,(ImgSize2+31):end);
        expo_mask_r2=expo_mask(56:120,(ImgSize2+31):end);
        expo_mask_r3=expo_mask(121:end,(ImgSize2+31):end);
        log_road_bright_r(ii)=mean(expo_mask_r1(expo_mask_r1~=0))+ ...
            mean(expo_mask_r2(expo_mask_r2~=0))+ ...
            mean(expo_mask_r3(expo_mask_r3~=0));
        log_road_bright_r(ii)=log_road_bright_r(ii)/3;
        log_hump(ii)=0; log_zebra(ii)=0;
        temp_p=0.3; temp_H=1.66; temp_td=230;
        %     log_xyimg{ii}=[]; % delet
        template_matching;
        win_r=16; win_c=30;
        conv_win=ones(win_r,win_c);
        ratio=(((1:(win_r/2))+win_r/2)/win_r)';
        ratio_top=repmat(ratio,1,150/2);
        ratio=((((win_r/2):-1:1)+win_r/2)/win_r)';
        ratio_btm=repmat(ratio,1,150/2);
% % % % %         tempimg=rightimg(1:2:end,:);
% % % % %         [~,A]=C_ipm_no_rotate2(tempimg,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);
% % % % %         Bori=A(1:5:end,1:2:end);
        pca_t=tic;
        tempimg=expo_mask(:,671:end);
        Aori=A;
        [~,A]=C_ipm_no_rotate2(tempimg,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);
        B=A(1:5:end,1:2:end);
        if ~isempty(ximg)
            ximg1_arrow=xyAsamp(1,1:4); yimg1_arrow=xyAsamp(2,1:4);
            min_x_arrow=max([1, min(ximg1_arrow)]);   max_x_arrow=min([size(B,2), max(ximg1_arrow)]);
            min_y_arrow=max([1, min(yimg1_arrow)]);             max_y_arrow=min([ImgSize1/2, max(yimg1_arrow)]);
            Bimgtmp=zeros(size(B,1),size(B,2));
            Bimgtmp(min_y_arrow:max_y_arrow,min_x_arrow:max_x_arrow)=B(min_y_arrow:max_y_arrow,min_x_arrow:max_x_arrow);
            [Bimgtmp_y, Bimgtmp_x]=find(Bimgtmp==1);
            if ~isempty(Bimgtmp_x)
                Bimgtmp_y=Bimgtmp_y-1;
                Bimgtmp_x=Bimgtmp_x-1;
                ximg1_arrow=ximg1_arrow-1;            yimg1_arrow=yimg1_arrow-1;
                C_in_poly_new(B, Bimgtmp_x, Bimgtmp_y, ximg1_arrow, yimg1_arrow);
            end
        end
        
        B(B==-1)=0; B_tip=B;
        B=B(1:4:end,1:2:end);
        B1=conv2(B,conv_win,'same');
        B1(1:(win_r/2),:)=B1(1:(win_r/2),:)./ratio_top;
        B1((end-floor(win_r/2)+1):end,:)=B1((end-floor(win_r/2)+1):end,:)./ratio_btm;
        B=[zeros(ceil(win_r/3),size(B,2)); B; zeros(ceil(win_r/3),size(B,2))];
        B1=[zeros(ceil(win_r/3),size(B,2)); B1;zeros(ceil(win_r/3),size(B,2))];
        [Btmp, Btmp2]=C_PCA_temp3_tmp(B,[win_r win_c],u(:,1:10),g(1:10,:),img_pca_mean,B1);
        % % %             [Btmp_line, Btmp2_line]=C_PCA_temp3_tmp(B,[win_r win_c],u_line(:,1:10),g_line(1:10,:),img_pca_mean_line,B1);
        
        log_minB(ii)=min(Btmp(:));
        [rr, cc]=find(Btmp<5.5);
        log_p(ii)=length(rr);
        [rr_min, cc_min]=find(Btmp==log_minB(ii));
        log_orien(ii)=Btmp2(rr_min(1), cc_min(1));
        %             Bori=[-ones(4*ceil(win_r/3),size(Bori,2)); Bori; -ones(4*ceil(win_r/3),size(Bori,2))];
        rr_tip=[];
        if ~isempty(rr)
            new_rr=rr-ceil(win_r/3);
            B_rr=4*new_rr-3;
            B_cc=2*cc-1;
            new_rr_min=rr_min-ceil(win_r/3);
            B_rr_min=4*(new_rr_min(1)-1)-3;
            B_cc_min=2*(cc_min(1)-1)-1;
            
            % % %                 B_rr=B_rr+4*ceil(win_r/3);
            % % %                 B_rr_min=B_rr_min+4*ceil(win_r/3);
            % % %                 indtmp=sub2ind(size(Btmp),rr,cc);
            % % %                 orien=Btmp2(indtmp)+1;
            % % %                 temp=img_pca_hi(:,orien);
            % % %                 [B_rr_min,B_cc_min,logdiff(ii)]=C_arrow_match(B_rr,B_cc,Bori,temp,[63 59]);
        elseif (isempty(rr) && log_p(ii-1)>0) || (isempty(rr) && log_p_tip(ii-1)>0)
            win_r_tip=27; win_c_tip=29;
            conv_win_tip=ones(win_r_tip,win_c_tip);
            B_tip_row_sum=sum(B_tip,2);
            row_ind=find(B_tip_row_sum~=0,1,'last');
            if ~isempty(row_ind)
                B_tip(row_ind:end,:)=repmat(B_tip(row_ind,:),size(B_tip,1)-row_ind+1,1);
            end
            B_tip=[zeros(ceil(win_r_tip/3),size(B_tip,2)); B_tip; zeros(ceil(win_r_tip/3),size(B_tip,2))];
            B1_tip=conv2(B_tip,conv_win_tip,'same');
            
%             B1_tip=[zeros(ceil(win_r_tip/3),size(B_tip,2)); B1_tip;zeros(ceil(win_r_tip/3),size(B_tip,2))];
            [Btmp_tip, Btmp2_tip]=C_PCA_temp3_tmp_tip(B_tip,[win_r_tip win_c_tip],u_tip(:,1:7),g_tip(1:7,:),img_pca_mean_tip,B1_tip);
            log_minB_tip(ii)=min(Btmp_tip(:));
            [rr_tip, cc_tip]=find(Btmp_tip<3.5);
            log_p_tip(ii)=length(rr_tip);
            [rr_min_tip, cc_min_tip]=find(Btmp_tip==log_minB_tip(ii));
            log_orien_tip(ii)=Btmp2_tip(rr_min_tip(1), cc_min_tip(1));
        else
            log_p_tip(ii)=0;
        end
        log_pca_t(ii)=toc(pca_t);
        % %             tic;
        % %             [cimtmp,r,c]=harris(B,1,2.0,3,0);
        % %             logt(ii)=toc;
        rimgtmp3d=imread(strcat(rightimg_dir,int2str(ii),'.jpg'));
        rimgtmp1d=double(rimgtmp3d(:,:,1));
        rimgtmp1d=Retify(rimgtmp1d,a1_right,a2_right,a3_right,a4_right,ind_1_right,ind_2_right,ind_3_right,ind_4_right,ind_new_right);
        rimgtmp3d(:,:,1)=rimgtmp1d;
        rimgtmp1d=double(rimgtmp3d(:,:,2));
        rimgtmp1d=Retify(rimgtmp1d,a1_right,a2_right,a3_right,a4_right,ind_1_right,ind_2_right,ind_3_right,ind_4_right,ind_new_right);
        rimgtmp3d(:,:,2)=rimgtmp1d;
        rimgtmp1d=double(rimgtmp3d(:,:,3));
        rimgtmp1d=Retify(rimgtmp1d,a1_right,a2_right,a3_right,a4_right,ind_1_right,ind_2_right,ind_3_right,ind_4_right,ind_new_right);
        rimgtmp3d(:,:,3)=rimgtmp1d;
        rimgtmp3d=rimgtmp3d(1:2:ImgSize1,:,:);
        imshow(rimgtmp3d,[]);
        if ~isempty(rr)
            vetex_rc=log_vetex{log_orien(ii)+1};
            B_vetex_r=4*(vetex_rc(1,:)+rr_min(1)-win_r/2-ceil(win_r/3))-3;
            B_vetex_c=2*(vetex_rc(2,:)+cc_min(1)-win_c/2)-1;
            A_vetex_r=5*B_vetex_r-4;
            A_vetex_c=2*B_vetex_c-1;
            xW_arrow=(150-A_vetex_c)/50; zW_arrow=(500-A_vetex_r)/50;
            xC_arrow=focal*xW_arrow./(temp_H*sin(temp_p)+cos(temp_p)*zW_arrow);
            yC_arrow=(zW_arrow*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_arrow*cos(temp_p));
            ximg_arrow=-xC_arrow+cc_x_right;
            yimg_arrow=-yC_arrow+cc_y;
            hold on;plot(ximg_arrow,yimg_arrow/2,'r');
            hold off;
        elseif isempty(rr) && ~isempty(rr_tip)
            vetex_rc_tip=log_vetex_tip{log_orien_tip(ii)+1};
            B_vetex_r_tip=vetex_rc_tip(1,:)+rr_min_tip(1)-win_r_tip/2-ceil(win_r_tip/3);
            B_vetex_c_tip=vetex_rc_tip(2,:)+cc_min_tip(1)-win_c_tip/2;
            A_vetex_r_tip=5*B_vetex_r_tip-4;
            A_vetex_c_tip=2*B_vetex_c_tip-1;
            xW_arrow_tip=(150-A_vetex_c_tip)/50; zW_arrow_tip=(500-A_vetex_r_tip)/50; zW_arrow_tip(zW_arrow_tip<0)=0;
            xC_arrow_tip=focal*xW_arrow_tip./(temp_H*sin(temp_p)+cos(temp_p)*zW_arrow_tip);
            yC_arrow_tip=(zW_arrow_tip*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_arrow_tip*cos(temp_p));
            ximg_arrow_tip=-xC_arrow_tip+cc_x_right;
            yimg_arrow_tip=-yC_arrow_tip+cc_y;
            hold on;plot(ximg_arrow_tip,yimg_arrow_tip/2,'b');
            hold off;
        end
        if ~isempty(ximg)
            hold on;plot(ximg,yimg/2,'g');
            hold off;
        end
        
% % % % %         imshow(Bori,[]);
% % % % %         if ~isempty(rr)
% % % % %             hold on; plot(B_cc,B_rr,'r.');
% % % % %             plot(B_cc_min(1),B_rr_min(1),'g.');
% % % % %             hold off;
% % % % %         end
% % % % %         if ~isempty(ximg)
% % % % %             hold on;plot([xyAsamp(1,:) xyAsamp(1,1)],[xyAsamp(2,:) xyAsamp(2,1)],'r');
% % % % %             hold off;
% % % % %         end
% % % % %         if ~isempty(rr_tip)
% % % % %             hold on; plot(cc_tip,rr_tip-ceil(win_r_tip/3),'b.');
% % % % %             plot(cc_min_tip(1),rr_min_tip(1)-ceil(win_r_tip/3),'c*');
% % % % %             hold off;
% % % % %         end
        drawnow;
        frame=getframe(gcf);
        im=frame2im(frame);
        imwrite(im,strcat(save_dir,'\3DArrow',int2str(ii),'.jpg'),'jpg');
        % % % % %             indtmp=find(Btmp==290);
        % % % % %             log290(ii)=length(indtmp);
        % % % % %             indtmp=find(Btmp==206);
        % % % % %             log206(ii)=length(indtmp);
        
        
    end
    P=find(log_minB<5.5);
    TP=intersect(P,True_ind);
    %     end
end