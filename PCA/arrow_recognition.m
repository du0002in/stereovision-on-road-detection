% if ~isempty(log_expos_box{ii-1})
%     expos_box_arrow=log_expos_box{ii-1};
% end
expos_box_arrow=expos_box;
expo_mask_arrow=C_in_region(img, expos_box_arrow{1}+40, expos_box_arrow{2}-40, expos_box_arrow{3}+ImgSize2+30+40, expos_box_arrow{4}+ImgSize2+30-40, expos_box_arrow{5}, expos_box_arrow{6});
expo_mask_arrow(1:25,:)=0;
tempimg=expo_mask_arrow(:,671:end);

if log_zebra(ii)==1
    ximg1=ximg_zebra(1:4); yimg1=yimg_zebra(1:4)/2;
    min_x=max([1, min(ximg1)]);   max_x=min([ImgSize2, max(ximg1)]);
    min_y=max([1, min(yimg1)]);   max_y=min([ImgSize1/2, max(yimg1)]);
    maxconvimgtmp=zeros(ImgSize1/2,ImgSize2);
    maxconvimgtmp(min_y:max_y,min_x:max_x)=tempimg(min_y:max_y,min_x:max_x);
    [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp~=0);
    if ~isempty(maxconvimgtmp_x)
        maxconvimgtmp_y=maxconvimgtmp_y-1;
        maxconvimgtmp_x=maxconvimgtmp_x-1;
        ximg1=ximg1-1;            yimg1=yimg1-1;
        C_in_poly_new(tempimg, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
    end
end

if log_hump(ii)==1
    ximg1=ximg_hump(1:4); yimg1=yimg_hump(1:4)/2;
    min_x=max([1, min(ximg1)]);   max_x=min([ImgSize2, max(ximg1)]);
    min_y=max([1, min(yimg1)]);   max_y=min([ImgSize1/2, max(yimg1)]);
    maxconvimgtmp=zeros(ImgSize1/2,ImgSize2);
    maxconvimgtmp(min_y:max_y,min_x:max_x)=tempimg(min_y:max_y,min_x:max_x);
    [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp~=0);
    if ~isempty(maxconvimgtmp_x)
        maxconvimgtmp_y=maxconvimgtmp_y-1;
        maxconvimgtmp_x=maxconvimgtmp_x-1;
        ximg1=ximg1-1;            yimg1=yimg1-1;
        C_in_poly_new(tempimg, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
    end
end

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
[rr_arrow, cc_arrow]=find(Btmp<5.5);
log_p(ii)=length(rr_arrow);
[rr_min, cc_min]=find(Btmp==log_minB(ii));
log_orien(ii)=Btmp2(rr_min(1), cc_min(1));
%             Bori=[-ones(4*ceil(win_r/3),size(Bori,2)); Bori; -ones(4*ceil(win_r/3),size(Bori,2))];
rr_tip=[];
if ~isempty(rr_arrow)
    new_rr=rr_arrow-ceil(win_r/3);
    B_rr=4*new_rr-3;
    B_cc=2*cc_arrow-1;
    new_rr_min=rr_min-ceil(win_r/3);
    B_rr_min=4*(new_rr_min(1)-1)-3;
    B_cc_min=2*(cc_min(1)-1)-1;
    
    % % %                 B_rr=B_rr+4*ceil(win_r/3);
    % % %                 B_rr_min=B_rr_min+4*ceil(win_r/3);
    % % %                 indtmp=sub2ind(size(Btmp),rr,cc);
    % % %                 orien=Btmp2(indtmp)+1;
    % % %                 temp=img_pca_hi(:,orien);
    % % %                 [B_rr_min,B_cc_min,logdiff(ii)]=C_arrow_match(B_rr,B_cc,Bori,temp,[63 59]);
elseif (isempty(rr_arrow) && log_p(ii-1)>0) || (isempty(rr_arrow) && log_p_tip(ii-1)>0)
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

if ~isempty(rr_arrow)
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
elseif isempty(rr_arrow) && ~isempty(rr_tip)
    vetex_rc_tip=log_vetex_tip{log_orien_tip(ii)+1};
    B_vetex_r_tip=vetex_rc_tip(1,:)+rr_min_tip(1)-win_r_tip/2-ceil(win_r_tip/3);
    B_vetex_c_tip=vetex_rc_tip(2,:)+cc_min_tip(1)-win_c_tip/2;
    A_vetex_r_tip=5*B_vetex_r_tip-4;
    A_vetex_c_tip=2*B_vetex_c_tip-1;
    xW_arrow_tip=(150-A_vetex_c_tip)/50; zW_arrow_tip=(500-A_vetex_r_tip)/50; 
    if sum(zW_arrow_tip<0)~=0
        if (zW_arrow_tip(1)<=0 && zW_arrow_tip(2)<=0)
            zW_arrow_tip(1)=0; zW_arrow_tip(2)=0;
        elseif (zW_arrow_tip(1)>0 && zW_arrow_tip(2)<=0)
            xW_arrow_tip(2)=xW_arrow_tip(1)+(xW_arrow_tip(1)-xW_arrow_tip(2))*zW_arrow_tip(1)/(zW_arrow_tip(2)-zW_arrow_tip(1));
            zW_arrow_tip(2)=0;
        end
        if (zW_arrow_tip(4)<=0 && zW_arrow_tip(3)<=0)
            zW_arrow_tip(4)=0; zW_arrow_tip(3)=0;
        elseif (zW_arrow_tip(4)>0 && zW_arrow_tip(3)<=0)
            xW_arrow_tip(3)=xW_arrow_tip(4)+(xW_arrow_tip(4)-xW_arrow_tip(3))*zW_arrow_tip(4)/(zW_arrow_tip(3)-zW_arrow_tip(4));
            zW_arrow_tip(3)=0;
        end
    end
    xC_arrow_tip=focal*xW_arrow_tip./(temp_H*sin(temp_p)+cos(temp_p)*zW_arrow_tip);
    yC_arrow_tip=(zW_arrow_tip*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_arrow_tip*cos(temp_p));
    ximg_arrow_tip=-xC_arrow_tip+cc_x_right;
    yimg_arrow_tip=-yC_arrow_tip+cc_y;
    ximg_arrow_tip_pre=ximg_arrow_tip;
    yimg_arrow_tip_pre=yimg_arrow_tip;
elseif isempty(rr_arrow) && isempty(rr_tip) && (log_p_tip(ii-1)>0 || log_p_tip(ii-2)>0)
% % % % %     letter_y_down_l_arr=(yimg_arrow_tip_pre(1)+yimg_arrow_tip_pre(2))*0.5;
% % % % %     letter_y_down_r_arr=(yimg_arrow_tip_pre(4)+yimg_arrow_tip_pre(3))*0.5;
% % % % %     letter_x_l_arr=(ximg_arrow_tip_pre(1)+ximg_arrow_tip_pre(2))*0.5;
% % % % %     letter_x_r_arr=(ximg_arrow_tip_pre(4)+ximg_arrow_tip_pre(3))*0.5;
% % % % %     esti_ximg_arr(1)=letter_x_l_arr;
% % % % %     esti_ximg_arr(2)=letter_x_r_arr;
% % % % %     esti_yimg_arr(1)=letter_y_down_l_arr;
% % % % %     esti_yimg_arr(2)=letter_y_down_r_arr;
% % % % %     esti_yimg_arr(3)=max([yimg_arrow_tip_pre(3)+letter_y_down_l_arr-yimg_arrow_tip_pre(1) 480]);
% % % % %     esti_yimg_arr(4)=max([yimg_arrow_tip_pre(4)+letter_y_down_r_arr-yimg_arrow_tip_pre(2) 480]);
% % % % %     esti_yimg_arr(5)=esti_yimg_arr(1);
% % % % %     esti_ximg_arr(3)=(esti_yimg_arr(3)-yimg_arrow_tip_pre(3))/(yimg_arrow_tip_pre(3)-yimg_arrow_tip_pre(4)) ...
% % % % %         *(ximg_arrow_tip_pre(3)-ximg_arrow_tip_pre(4))+ximg_arrow_tip_pre(3);
% % % % %     esti_ximg_arr(4)=(esti_yimg_arr(4)-yimg_arrow_tip_pre(2))/(yimg_arrow_tip_pre(2)-yimg_arrow_tip_pre(1)) ...
% % % % %         *(ximg_arrow_tip_pre(2)-ximg_arrow_tip_pre(1))+ximg_arrow_tip_pre(2);
% % % % %     esti_ximg_arr(5)=esti_ximg_arr(1);
% % %     esti_ximg_arr=ximg_arrow_tip_pre;
% % %     esti_yimg_arr=yimg_arrow_tip_pre;
% % %     esti_yimg_arr([1 4])=(yimg_arrow_tip_pre([1, 4])+ImgSize1)*0.5;
% % %     esti_yimg_arr([2 3])=yimg_arrow_tip_pre([2, 3])+(esti_yimg_arr([1 4])-yimg_arrow_tip_pre([1, 4]));
% % %     esti_yimg_arr(5)=esti_yimg_arr(1); 
% % %     ximg_arrow_tip_est=esti_ximg_arr; yimg_arrow_tip_est=esti_yimg_arr;
% % % %     dist_move=sum(loop_v.*loop_t);
% % % %     zW_arrow_tip=zW_arrow_tip-dist_move;
% % % %     zW_arrow_tip(zW_arrow_tip<0)=0;
% % % %     xC_arrow_tip=focal*xW_arrow_tip./(temp_H*sin(temp_p)+cos(temp_p)*zW_arrow_tip);
% % % %     yC_arrow_tip=(zW_arrow_tip*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_arrow_tip*cos(temp_p));
% % % %     ximg_arrow_tip_est=-xC_arrow_tip+cc_x_right;
% % % %     yimg_arrow_tip_est=-yC_arrow_tip+cc_y;
    
    new_states=C_Predict_P([0;0;pi/2],loop_t,loop_phi,loop_v);
    d_theta_hump=-new_states(3)+pi/2;
    R_hump=[cos(d_theta_hump) -sin(d_theta_hump); sin(d_theta_hump) cos(d_theta_hump)];
    T_hump=-R_hump*[new_states(1); new_states(2)];
    RT_hump=[R_hump T_hump; 0 0 1];
    xzW_arrow_tip=RT_hump*[xW_arrow_tip;zW_arrow_tip+L_wheel_cam;ones(1,length(xW_arrow_tip))];
    xW_arrow_tip=xzW_arrow_tip(1,:); zW_arrow_tip=xzW_arrow_tip(2,:); 
    zW_arrow_tip=zW_arrow_tip-L_wheel_cam;
    if sum(zW_arrow_tip<0)~=0
        if (zW_arrow_tip(1)<=0 && zW_arrow_tip(2)<=0)
            zW_arrow_tip(1)=0; zW_arrow_tip(2)=0;
        elseif (zW_arrow_tip(1)>0 && zW_arrow_tip(2)<=0)
            xW_arrow_tip(2)=xW_arrow_tip(1)+(xW_arrow_tip(1)-xW_arrow_tip(2))*zW_arrow_tip(1)/(zW_arrow_tip(2)-zW_arrow_tip(1));
            zW_arrow_tip(2)=0;
        end
        if (zW_arrow_tip(4)<=0 && zW_arrow_tip(3)<=0)
            zW_arrow_tip(4)=0; zW_arrow_tip(3)=0;
        elseif (zW_arrow_tip(4)>0 && zW_arrow_tip(3)<=0)
            xW_arrow_tip(3)=xW_arrow_tip(4)+(xW_arrow_tip(4)-xW_arrow_tip(3))*zW_arrow_tip(4)/(zW_arrow_tip(3)-zW_arrow_tip(4));
            zW_arrow_tip(3)=0;
        end
    end
    xC_arrow_tip=focal*xW_arrow_tip./(temp_H*sin(temp_p)+cos(temp_p)*zW_arrow_tip);
    yC_arrow_tip=(zW_arrow_tip*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_arrow_tip*cos(temp_p));
    ximg_arrow_tip_est=-xC_arrow_tip+cc_x_right;
    yimg_arrow_tip_est=-yC_arrow_tip+cc_y;
end