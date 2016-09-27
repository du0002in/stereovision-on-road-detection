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
[ipm_img,~]=C_ipm_no_rotate2(tempimg,temp_p,temp_H,focal,temp_td,cc_x_right,cc_y);
log_hump(ii)=0;
B_hump=ipm_img(1:2:end,1:2:end);
sum_B_hump_row=sum(B_hump,2);
convB_hump=conv2(B_hump,ones(1,5)/5,'same');
ridge_hump=C_new_Ridge_zebra(convB_hump,5,sum_B_hump_row);
sum_ridge_hump=sum(ridge_hump,2);
startrow_hump=1; endrow_hump=1;
current_hump_row=[];
for row_h=1:2:size(B_hump,1)
    if sum_ridge_hump(row_h)>=4
        b=find(ridge_hump(row_h,:)==1);
        lead0_ind=max([(b(1)-6) 1]);
        end0_ind=min([b(end)+6 size(B_hump,2)]);
        if (end0_ind-lead0_ind>=40)
            a=convB_hump(row_h,lead0_ind:end0_ind);
            a_ori=a;
            a=a-mean(a);
            locs=find(ridge_hump(row_h,lead0_ind:end0_ind)==1);
            if length(locs)>=4
                locs_diff=locs(2:end)-locs(1:(end-1));
                locs_diff_valid=(locs_diff>=10 & locs_diff<=18)*1.0;
                if (sum(locs_diff_valid)<3)
                    hump_status=0;
                else
                    one_ind=C_find_connect_ones(locs_diff_valid,3);
                    dumy_counter=one_ind(2)-one_ind(1)+1;
                    if dumy_counter>=3
                        hump_status=1;
                        lead0_ind=max([(locs(one_ind(1))-6) 1]);
                        end0_ind=min([locs(one_ind(2)+1)+6 size(B_hump,2)]);
                        a=a(lead0_ind:end0_ind);
                        a=a-mean(a);
                    else
                        hump_status=0;
                    end
                end
            else
                hump_status=0;
            end
            if hump_status==1
                a=fft(a,NFFT_hump)/length(a);
                a=a(1:(NFFT_hump/2+1));
                ma=abs(a);
                [ma_max,ma_ind]=max(ma);
                if (fre_seq_hump(ma_ind)>=9.5 && fre_seq_hump(ma_ind)<=16) && (ma_max>=7)
                    endrow_hump=row_h;
                    current_hump_row=[current_hump_row endrow_hump:(endrow_hump+1)];
                else
                    for j=3:dumy_counter
                        lead0_ind=max([(locs(one_ind(1)+j-3)-6) 1]);
                        end0_ind=min([locs(one_ind(1)+j)+6 size(B_hump,2)]);
                        a=a_ori(lead0_ind:end0_ind);
                        a=a-mean(a);
                        a=fft(a,NFFT_hump)/length(a);
                        a=a(1:(NFFT_hump/2+1));
                        ma=abs(a);
                        [ma_max,ma_ind]=max(ma);
                        if (fre_seq_hump(ma_ind)>=9.0 && fre_seq_hump(ma_ind)<=16) && (ma_max>=7)
                            endrow_hump=row_h;
                            current_hump_row=[current_hump_row endrow_hump:(endrow_hump+1)];
                            break;
                        end
                    end
                end
            end
        end
    end
    if (endrow_hump-startrow_hump)==2
        startrow_hump=endrow_hump;
    else
        startrow_hump=row_h;
    end
end
if ~isempty(current_hump_row)
    current_hump_row=C_group_rows(current_hump_row,5);
    current_hump_row(current_hump_row==0)=[];
    if ~isempty(current_hump_row)
        log_hump(ii)=1;
        if current_hump_row(end)<98
            current_hump_row(1)=current_hump_row(end)-98;
        end
    end
end
if log_hump(ii)==1
    B_vetex_hump_c=[1 size(B_hump,2) size(B_hump,2) 1 1];
    B_vetex_hump_r=[current_hump_row(1) current_hump_row(1) current_hump_row(end) current_hump_row(end) current_hump_row(1)];
    A_vetex_hump_c=2*B_vetex_hump_c-1;
    A_vetex_hump_r=2*B_vetex_hump_r-1;
    xW_hump=(150-A_vetex_hump_c)/50; zW_hump=(500-A_vetex_hump_r)/50;
    xC_hump=focal*xW_hump./(temp_H*sin(temp_p)+cos(temp_p)*zW_hump);
    yC_hump=(zW_hump*sin(temp_p)*focal-temp_H*focal*cos(temp_p))./(temp_H*sin(temp_p)+zW_hump*cos(temp_p));
    ximg_hump=-xC_hump+cc_x_right;
    yimg_hump=-yC_hump+cc_y;
end
