%template_matching
tempimg=rightimg(1:2:end,:);
if ~isempty(ximg)
    if ~isempty(log_xyimg{ii-1})
        ab=[ximg(1) 1; ximg(2) 1]\[yimg(1); yimg(2)];
        letter_ytmp=(ab(1)*(1:size(tempimg,2))+ab(2))/2-1;
        C_remove_pixel_above_line(tempimg,letter_ytmp);
    end
end

log_letter_status(ii)=0;
ximg=[]; yimg=[]; log_xyimg{ii}=[];

if (log_hump(ii-1)==1 || log_hump(ii-2)==1 || log_zebra(ii-1)==1 || log_zebra(ii-2)==1)
    return;
end

% [ipm_img,A]=C_ipm_no_rotate(tempimg,temp_p,temp_H,focal,xxc,yyc,temp_td);
[ipm_img,A]=C_ipm_no_rotate2(tempimg,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);
B=A(1:5:end,1:2:end);
Bp=B; Bp(B==-1)=0;
conv_B=C_conv_uni_accurate(Bp,25,53);
Bp2=B; Bp2(B==1)=0; Bp2=Bp2+1;
conv_Bp2=C_conv_uni_accurate(Bp2,25,53);
indtmp=(conv_Bp2>833 & conv_Bp2<1325);
conv_B(indtmp)=conv_B(indtmp)./conv_Bp2(indtmp)*1325;
C=conv_B;
% C(conv_B<=411 | conv_B>=795)=0;C(conv_B>411 & conv_B<795)=1;
C(conv_B<=411 | conv_B>=663)=0;C(conv_B>411 & conv_B<663)=1;
C(1:4,:)=0; C(end-3:end,:)=0;
ind=find(C==1)-1;
exp_C=C_expan_img(C,ind,25,53);
filtered_img=zeros(size(B,1),size(B,2),13);
for i=1:length(log_mask)
    filtered_img(:,:,i)=conv2(B,log_mask{i},'same')/sum(log_mask{i}(:));
end
[max_img, max_ind]=max(filtered_img,[],3);
indtmp=(max_img<0.75 | exp_C==0);
max_ind(indtmp)=0;
max_img(indtmp)=0;
for i=1:13
    ind=find(max_ind==i);
    log_ind{i}=ind;
    log_num(i)=length(ind);
end
[sort_log_num, seq]=sort(log_num,'descend');
inddumy=find(sort_log_num==0,1);
if inddumy==1
    return;
elseif isempty(inddumy)
    inddumy=length(sort_log_num);
end

% log_letter_dumy(ii)=1;
% log_ind_len(ii)=0;
for theta_ind=1:2
    temp_d_theta=theta(seq(theta_ind));
    
    [rotBp, move]=C_rot_img2(B,temp_d_theta);
    rotBp2=rotBp;
    rotBp2(rotBp~=3)=0; rotBp2(rotBp==3)=1;
    conv_B=C_conv_uni_accurate(rotBp2,25,53);
    rotBp3=rotBp;
    rotBp3(rotBp==1)=0;
    rotBp3(rotBp==2 | rotBp==3)=1;
    conv_Bp3=C_conv_uni_accurate(rotBp3,25,53);
    indtmp=(conv_Bp3>883 & conv_Bp3<1325);
    conv_B(indtmp)=conv_B(indtmp)./conv_Bp3(indtmp)*1325;
    C=conv_B;
    %     C(conv_B<=464 | conv_B>=795)=0;C(conv_B>464 & conv_B<795)=1;
    C(conv_B<=464 | conv_B>=663)=0;C(conv_B>464 & conv_B<663)=1;
    C(1:4,:)=0; C((end-3):end,:)=0;
    ind=find(C==1);
%     log_ind_len(ii)=max([log_ind_len(ii) length(ind)]);
    %     toc(ss)
    if isempty(ind)
        continue;
    end
    ind=ind(1:2:length(ind));
    rotBp=rotBp-2;
    counter=0;
    template_order;
    for i=itmp
        temp=roadmarkingtemplates{i};
        width=letter_width{i};
        cum_width=cumsum(width); cum_width=[0 cum_width];
        temp_hight=size(temp,1);
        if i==1
            log_atrribute='AHEAD';
        elseif i==4
            log_atrribute='HUMP';
        elseif i==5
            log_atrribute='SLOW';
        elseif i==6
            log_atrribute='STRIPS';
        elseif i==7
            log_atrribute='X-ING';
        elseif i==3
            log_atrribute='BUS';
        elseif i==10
            log_atrribute='SLOW';
        end
        if i==1
            jt=[3 4 1];
        elseif (i==10 || i==11)
            jt=2;
        else
            jt=[3 2];
        end
        for j=jt
            template=[];
            for k=1:length(width)
                template=[template temp(:,(cum_width(k)+1):cum_width(k+1)) -ones(temp_hight,j)];
            end
            template=double(template);
            template(:,(end-j+1):end)=[];
            counter=counter+1;
%             [log_min, log_min_ind]=min(C_L1_norm(rotBp,ind,template));
            [log_min, log_min_ind]=C_L1_norm2(rotBp,ind,template);
            if (i~=7 && log_min<0.30) || (i==7 && log_min<0.25)
                min_ind=ind(log_min_ind);
                [y, x]=ind2sub(size(rotBp2),min_ind);
                xrotBp2=[x-size(template,2)/2 x+size(template,2)/2 x+size(template,2)/2 x-size(template,2)/2]-1;
                yrotBp2=[y-size(template,1)/2 y-size(template,1)/2 y+size(template,1)/2 y+size(template,1)/2]-1;
                R=[cos(temp_d_theta) -sin(temp_d_theta)*2.5; sin(temp_d_theta)*0.4 cos(temp_d_theta)];
                xyAsamp=R\[xrotBp2-move(1); yrotBp2-move(2)];
                xyA=[xyAsamp(1,:)*2; xyAsamp(2,:)*5];
                xA=xyA(1,:)+1; yA=xyA(2,:)+1;
                xW_letter=(150-xA)/50; zW_letter=(500-yA)/50;
                xC_letter=focal*xW_letter./(temp_H*sin(temp_p)+cos(temp_p)*zW_letter);
                yC_letter=(zW_letter*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_letter*cos(temp_p));
                ximg=-xC_letter+cc_x_right;
                yimg=-yC_letter+cc_y;
                ximg=[ximg ximg(1)]; yimg=[yimg yimg(1)];
                total_counter=total_counter+1;
                log_data_name{total_counter}=log_atrribute;
                log_data_j(total_counter)=j;
                log_data_img(total_counter)=ii;
                log_data_mis(total_counter)=log_min;
                log_letter_status(ii)=i;
                log_xyimg{ii}=[ximg; yimg];
                break;
            end
        end
        if log_min<0.30
            break;
        end
    end
    if log_min<0.30
        break;
    end
end

% % if (isempty(ximg) && ~isempty(log_xyimg{ii-1}))
% % % % % % %     if max(log_xyimg{ii-1}(2,:))>300
% % % % % % %         letter_y_down_l=(log_xyimg{ii-1}(2,1)+log_xyimg{ii-1}(2,4))*0.5;
% % % % % % %         letter_y_down_r=(log_xyimg{ii-1}(2,2)+log_xyimg{ii-1}(2,3))*0.5;
% % % % % % %         letter_x_l=(log_xyimg{ii-1}(1,1)+log_xyimg{ii-1}(1,4))*0.5;
% % % % % % %         letter_x_r=(log_xyimg{ii-1}(1,2)+log_xyimg{ii-1}(1,3))*0.5;
% % % % % % %         esti_ximg(1)=letter_x_l;
% % % % % % %         esti_ximg(2)=letter_x_r;
% % % % % % %         esti_yimg(1)=letter_y_down_l;
% % % % % % %         esti_yimg(2)=letter_y_down_r;
% % % % % % %         esti_yimg(3)=max([log_xyimg{ii-1}(2,3)+letter_y_down_l-log_xyimg{ii-1}(2,1) 480]);
% % % % % % %         esti_yimg(4)=max([log_xyimg{ii-1}(2,4)+letter_y_down_r-log_xyimg{ii-1}(2,2) 480]);
% % % % % % %         esti_yimg(5)=esti_yimg(1);
% % % % % % %         esti_ximg(3)=(esti_yimg(3)-log_xyimg{ii-1}(2,3))/(log_xyimg{ii-1}(2,3)-log_xyimg{ii-1}(2,2)) ...
% % % % % % %             *(log_xyimg{ii-1}(1,3)-log_xyimg{ii-1}(1,2))+log_xyimg{ii-1}(1,3);
% % % % % % %         esti_ximg(4)=(esti_yimg(4)-log_xyimg{ii-1}(2,4))/(log_xyimg{ii-1}(2,4)-log_xyimg{ii-1}(2,1)) ...
% % % % % % %             *(log_xyimg{ii-1}(1,4)-log_xyimg{ii-1}(1,1))+log_xyimg{ii-1}(1,4);
% % % % % % %         esti_ximg(5)=esti_ximg(1);
% % % % % % %         ximg=esti_ximg; yimg=esti_yimg;
% % % % % % %     end
% %     new_states=C_Predict_P([0;0;pi/2],loop_t,loop_phi,loop_v);
% %     d_theta_letter=-new_states(3)+pi/2;
% %     R_letter=[cos(d_theta_letter) -sin(d_theta_letter); sin(d_theta_letter) cos(d_theta_letter)];
% %     T_letter=-R_letter*[new_states(1); new_states(2)];
% %     RT_letter=[R_letter T_letter; 0 0 1];
% %     xzW_letter=RT_letter*[xW_letter;zW_letter;ones(1,length(xW_letter))];
% %     xW_letter=xzW_letter(1,:); zW_letter=xzW_letter(2,:);
% %     xC_letter=focal*xW_letter./(temp_H*sin(temp_p)+cos(temp_p)*zW_letter);
% %     yC_letter=(zW_letter*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_letter*cos(temp_p));
% %     ximg=-xC_letter+cc_x_right;
% %     yimg=-yC_letter+cc_y;
% %     ximg=[ximg ximg(1)]; yimg=[yimg yimg(1)];
% % % % end