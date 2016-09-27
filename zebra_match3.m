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
% %for Image Collection(new Cam)8
log_ite_s=[3 284 461 3];
log_ite_e=[378 284 857 911];

%for Image Collection(new Cam)9
% log_ite_s=[101 703 461 3];
% log_ite_e=[530 1131 857 911];

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
    ImgSize1=480; ImgSize2=640; NFFT_zebra=2^nextpow2(ImgSize2);
    fre_seq_zebra=ImgSize2/2*linspace(0,1,NFFT_zebra/2+1);
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
    load 'zebra_gt.mat';
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
    %     KickOff_PF2;
    counter_zebra=2;
    log_abnorm=zeros(1,4);
    
    B_zebra=zeros(225,150);
    NFFT_zebra=2^nextpow2(size(B_zebra,2));
    fre_seq_zebra=size(B_zebra,2)/2*linspace(0,1,NFFT_zebra/2+1);
    for ii=[135:143, 284:287, 579:586, 733:740]
%     for ii=ite_s:ite_e
        counter_zebra=2;
        log_zebra(ii)=0;
        zebra_rows=[];
        if ii==138
            a=0;
        end
        in_start=tic;
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
        
        
        temp_p=0.3; temp_H=1.66; temp_td=230;
        tempimg=rightimg(1:2:end,:);
        [ipm_img,A]=C_ipm_no_rotate2(tempimg,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);
        
        z_t=tic;
        B_zebra=ipm_img(1:2:end,1:2:end); Bori=B_zebra;
        sum_B_zebra_row=sum(B_zebra,2);
        convB_zebra=conv2(B_zebra,ones(1,16)/16,'same');
        ridge_zebra_P=C_new_Ridge_zebra(convB_zebra,12,sum_B_zebra_row);
        sum_B_zebra_row=sum(ridge_zebra_P,2);
        ridge_zebra_N=C_new_Ridge_zebra(-convB_zebra,12,sum_B_zebra_row);
        ridge_zebra=ridge_zebra_P+ridge_zebra_N;
        ridge_zebra(ridge_zebra>1)=1;
        sum_ridge_zebra=sum(ridge_zebra,2);
        startrow_zebra=1; endrow_zebra=1;
        loop_t=tic;
        for row_zebra=1:2:size(B_zebra,1)
            if row_zebra==39 || row_zebra==56
                aa=0;
            end
            if sum_ridge_zebra(row_zebra)>=3
                b=find(ridge_zebra(row_zebra,:)==1);
                lead0_ind=max([(b(1)-8) 1]);
                end0_ind=min([b(end)+8 size(B_zebra,2)]);
                if (end0_ind-lead0_ind>=25)
                    a=convB_zebra(row_zebra,lead0_ind:end0_ind);
                    a=a-mean(a);
                    locs=find(ridge_zebra(row_zebra,lead0_ind:end0_ind)==1);
                    locs(abs(a(locs)-mean(a))<=10)=[];
                    if length(locs)>=2
                        locs_diff=locs(2:end)-locs(1:(end-1));
                        locs_diff_valid=(locs_diff>=11 & locs_diff<=20)*1.0;
                        if (sum(locs_diff_valid)<3 && row_zebra<=127) || (sum(locs_diff_valid)<2 && row_zebra>127)
                            zebra_status=0;
                        else
                            if row_zebra<=127
                                one_ind=C_find_connect_ones(locs_diff_valid,3);
                                dumy_counter=one_ind(2)-one_ind(1)+1;
                                if dumy_counter>=3
                                    zebra_status=1;
                                    lead0_ind=max([(locs(one_ind(1))-8) 1]);
                                    end0_ind=min([locs(one_ind(2)+1)+8 size(B_zebra,2)]);
                                    a=a(lead0_ind:end0_ind);
                                    a=a-mean(a);
                                else
                                    zebra_status=0;
                                end
                            else
                                one_ind=C_find_connect_ones(locs_diff_valid,2);
                                dumy_counter=one_ind(2)-one_ind(1)+1;
                                if dumy_counter>=2
                                    zebra_status=1;
                                    lead0_ind=max([(locs(one_ind(1))-8) 1]);
                                    end0_ind=min([locs(one_ind(2)+1)+8 size(B_zebra,2)]);
                                    a=a(lead0_ind:end0_ind);
                                    a=a-mean(a);
                                else
                                    zebra_status=0;
                                end
                            end
                            
                        end
                    else
                        zebra_status=0;
                    end
                    if zebra_status==1
                        a=fft(a,NFFT_zebra)/length(a);
                        a=a(1:(NFFT_zebra/2+1));
                        ma=abs(a);
                        
                        [ma_max,ma_ind]=max(ma);
                        if (fre_seq_zebra(ma_ind)>=4 && fre_seq_zebra(ma_ind)<=6.5) && (ma_max>=7)
                            endrow_zebra=row_zebra;
                        end
                    end
                end
            end
            if (endrow_zebra-startrow_zebra)==2
                counter_zebra=counter_zebra+2;
                startrow_zebra=endrow_zebra;
            else
                if counter_zebra>=10
                    B_zebra((endrow_zebra-counter_zebra+1):endrow_zebra,:)=255;
                    log_zebra(ii)=1;
                    zebra_rows=[zebra_rows;(endrow_zebra-counter_zebra+1) endrow_zebra];
                end
                counter_zebra=1;
                startrow_zebra=row_zebra;
            end
        end
        log_loopt(ii)=toc(loop_t);
        ridge2_zebra=zeros(size(ridge_zebra));
        if log_zebra(ii)==1
            for i=1:size(zebra_rows,1)
                ridge2_zebra(zebra_rows(i,1):zebra_rows(i,2),:)=ridge_zebra(zebra_rows(i,1):zebra_rows(i,2),:);
            end
            [rr_zebra,cc_zebra]=find(ridge2_zebra==1);
            [line_in, line_para_zebra]=C_fit_single_line(rr_zebra,cc_zebra);
            a_l_zebra=line_para_zebra(1); b_l_zebra=line_para_zebra(2); c_l_zebra=line_para_zebra(3);
            if (a_l_zebra*b_l_zebra>0)
                c_l_left=-a_l_zebra*cc_zebra(1)-b_l_zebra*rr_zebra(1);
                c_l_right=-a_l_zebra*cc_zebra(end)-b_l_zebra*rr_zebra(end);
                a_l1=a_l_zebra; b_l1=b_l_zebra; c_l1=c_l_left+8.0*sqrt(a_l_zebra^2+b_l_zebra^2);
                a_l3=a_l_zebra; b_l3=b_l_zebra; c_l3=c_l_right-8.0*sqrt(a_l_zebra^2+b_l_zebra^2);
                if a_l_zebra~=0
                    x_left_most=(-c_l1-b_l1*zebra_rows(end,2))/a_l1;
                    x_right_most=(-c_l3-b_l3*zebra_rows(1,1))/a_l3;
                else
                    x_left_most=cc_zebra(1)-8;
                    x_right_most=cc_zebra(end)+8;
                end
                if a_l_zebra~=0
                    b_l4=1.0;
                    a_l4=-b_l1*b_l4/a_l1;
                    c_l4=-(a_l4*x_left_most+b_l4*zebra_rows(end,2));
                    
                    b_l2=1.0;
                    a_l2=-b_l3*b_l2/a_l3;
                    c_l2=-(a_l2*x_right_most+b_l2*zebra_rows(1,1));
                else
                    b_l4=0; a_l4=1;
                    c_l4=-a_l4*x_left_most;
                    
                    b_l2=0; a_l2=1;
                    c_l2=-a_l2*x_right_most;
                end
                x41=x_left_most;    y41=zebra_rows(end,2);
                x23=x_right_most;   y23=zebra_rows(1,1);
                [x12, y12]=line_intersect([a_l1 b_l1 c_l1],[a_l2 b_l2 c_l2]);
                [x34, y34]=line_intersect([a_l3 b_l3 c_l3],[a_l4 b_l4 c_l4]);
            else
                c_l_left=-a_l_zebra*cc_zebra(1)-b_l_zebra*rr_zebra(1);
                c_l_right=-a_l_zebra*cc_zebra(end)-b_l_zebra*rr_zebra(end);
                a_l1=a_l_zebra; b_l1=b_l_zebra; c_l1=c_l_left+8.0*sqrt(a_l_zebra^2+b_l_zebra^2);
                a_l3=a_l_zebra; b_l3=b_l_zebra; c_l3=c_l_right-8.0*sqrt(a_l_zebra^2+b_l_zebra^2);
                if a_l_zebra~=0
                    x_left_most=(-c_l1-b_l_zebra*zebra_rows(1,1))/a_l_zebra;
                    x_right_most=(-c_l3-b_l_zebra*zebra_rows(end,2))/a_l_zebra;
                else
                    x_left_most=cc_zebra(1)-8;
                    x_right_most=cc_zebra(end)+8;
                end
                if a_l_zebra~=0
                    b_l4=1.0;
                    a_l4=-b_l1*b_l4/a_l1;
                    c_l4=-(a_l4*x_right_most+b_l4*zebra_rows(end,2));
                    
                    b_l2=1.0;
                    a_l2=-b_l3*b_l2/a_l3;
                    c_l2=-(a_l2*x_left_most+b_l2*zebra_rows(1,1));
                else
                    b_l4=0; a_l4=1;
                    c_l4=-a_l4*x_right_most;
                    
                    b_l2=0; a_l2=1;
                    c_l2=-a_l2*x_left_most;
                end
                x12=x_left_most;    y12=zebra_rows(1,1);
                x34=x_right_most;   y34=zebra_rows(end,2);
                [x23, y23]=line_intersect([a_l2 b_l2 c_l2],[a_l3 b_l3 c_l3]);
                [x41, y41]=line_intersect([a_l1 b_l1 c_l1],[a_l4 b_l4 c_l4]);
            end
        end
        log_z_t(ii)=toc(z_t);
        
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
        
        %         imshow([Bori,B,ridge2*255],[]);
        %         hold on;plot(cc(line_in==1)+300,rr(line_in==1),'r.'); hold off;
        % %         if log_zebra(ii)==1
        % %             hold on;plot([x12 x23 x34 x41 x12],[y12 y23 y34 y41 y12],'r');
        % %         end
        if log_zebra(ii)==1
            B_vetex_zebra_c=[x12 x23 x34 x41 x12];
            B_vetex_zebra_r=[y12 y23 y34 y41 y12];
            A_vetex_zebra_c=2*B_vetex_zebra_c-1;
            A_vetex_zebra_r=2*B_vetex_zebra_r-1;
            xW_zebra=(150-A_vetex_zebra_c)/50; zW_zebra=(500-A_vetex_zebra_r)/50;
            xC_zebra=focal*xW_zebra./(temp_H*sin(temp_p)+cos(temp_p)*zW_zebra);
            yC_zebra=(zW_zebra*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_zebra*cos(temp_p));
            ximg_zebra=-xC_zebra+cc_x_right;
            yimg_zebra=-yC_zebra+cc_y;
            hold on;plot(ximg_zebra, yimg_zebra/2, 'r', 'LineWidth', 2); hold off;
        end
        drawnow;
        frame=getframe(gcf);
        im=frame2im(frame);
        imwrite(im,strcat(save_dir,'\3DZebra',int2str(ii),'.jpg'),'jpg');
        % %         log_TP_counter(ii)=TP_counter;
        % %         log_FP_counter(ii)=FP_counter;
    end
    %     end
end