tempimg=rightimg(1:2:end,:);
if ~isempty(ximg)
    ximg1=ximg(1:4); yimg1=yimg(1:4)/2;
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
% % % % % 
% % % % % if ~isempty(rr_arrow)
% % % % %     ximg1=ximg_arrow(1:4); yimg1=yimg_arrow(1:4)/2;
% % % % %     min_x=max([1, min(ximg1)]);   max_x=min([ImgSize2, max(ximg1)]);
% % % % %     min_y=max([1, min(yimg1)]);   max_y=min([ImgSize1/2, max(yimg1)]);
% % % % %     maxconvimgtmp=zeros(ImgSize1/2,ImgSize2);
% % % % %     maxconvimgtmp(min_y:max_y,min_x:max_x)=tempimg(min_y:max_y,min_x:max_x);
% % % % %     [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp~=0);
% % % % %     if ~isempty(maxconvimgtmp_x)
% % % % %         maxconvimgtmp_y=maxconvimgtmp_y-1;
% % % % %         maxconvimgtmp_x=maxconvimgtmp_x-1;
% % % % %         ximg1=ximg1-1;            yimg1=yimg1-1;
% % % % %         C_in_poly_new(tempimg, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
% % % % %     end
% % % % % elseif isempty(rr_arrow) && ~isempty(rr_tip)
% % % % %     ximg1=ximg_arrow_tip(1:4); yimg1=yimg_arrow_tip(1:4)/2;
% % % % %     min_x=max([1, min(ximg1)]);   max_x=min([ImgSize2, max(ximg1)]);
% % % % %     min_y=max([1, min(yimg1)]);   max_y=min([ImgSize1/2, max(yimg1)]);
% % % % %     maxconvimgtmp=zeros(ImgSize1/2,ImgSize2);
% % % % %     maxconvimgtmp(min_y:max_y,min_x:max_x)=tempimg(min_y:max_y,min_x:max_x);
% % % % %     [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp~=0);
% % % % %     if ~isempty(maxconvimgtmp_x)
% % % % %         maxconvimgtmp_y=maxconvimgtmp_y-1;
% % % % %         maxconvimgtmp_x=maxconvimgtmp_x-1;
% % % % %         ximg1=ximg1-1;            yimg1=yimg1-1;
% % % % %         C_in_poly_new(tempimg, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
% % % % %     end
% % % % % end

counter_zebra=2;
log_zebra(ii)=0;
zebra_rows=[];
[ipm_img,A]=C_ipm_no_rotate2(tempimg,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);
B_zebra=ipm_img(1:2:end,1:2:end);
sum_B_zebra_row=sum(B_zebra,2);
convB_zebra=conv2(B_zebra,ones(1,16)/16,'same');
ridge_zebra_P=C_new_Ridge_zebra(convB_zebra,12,sum_B_zebra_row);
sum_B_zebra_row=sum(ridge_zebra_P,2);
ridge_zebra_N=C_new_Ridge_zebra(-convB_zebra,12,sum_B_zebra_row);
ridge_zebra=ridge_zebra_P+ridge_zebra_N;
ridge_zebra(ridge_zebra>1)=1;
sum_ridge_zebra=sum(ridge_zebra,2);
startrow_zebra=1; endrow_zebra=1;
for row_zebra=1:2:size(B_zebra,1)
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
if log_zebra(ii)==1
    ridge2=zeros(size(ridge_zebra));
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
    B_vetex_zebra_c=[x12 x23 x34 x41 x12];
    B_vetex_zebra_r=[y12 y23 y34 y41 y12];
    A_vetex_zebra_c=2*B_vetex_zebra_c-1;
    A_vetex_zebra_r=2*B_vetex_zebra_r-1;
    xW_zebra=(150-A_vetex_zebra_c)/50; zW_zebra=(500-A_vetex_zebra_r)/50;
    xC_zebra=focal*xW_zebra./(temp_H*sin(temp_p)+cos(temp_p)*zW_zebra);
    yC_zebra=(zW_zebra*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_zebra*cos(temp_p));
    ximg_zebra=-xC_zebra+cc_x_right;
    yimg_zebra=-yC_zebra+cc_y;
end




