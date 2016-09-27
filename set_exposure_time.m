% src1.Exposure=exposure_time_l;
% src2.Exposure=exposure_time_r;
% trigger([vid1 vid2]);
% leftimg=getdata(vid1);
% rightimg=getdata(vid2);
leftimg_tmp=leftimg; rightimg_tmp=rightimg;
leftimg=double(max(leftimg_tmp,[],3));
rightimg=double(max(rightimg_tmp,[],3));
img=[leftimg(1:2:ImgSize1,:) zeros(ImgSize1/2, 30) rightimg(1:2:ImgSize1,:)];
expo_mask=C_in_region(img, expos_box{1}, expos_box{2}, expos_box{3}+ImgSize2+30, expos_box{4}+ImgSize2+30, expos_box{5}, expos_box{6});
% expo_mask(1:25,:)=0;
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
if log_road_bright_l(ii)>150 || log_road_bright_l(ii)<100
    exposure_time_l_tmp=exposure_time_l*133/log_road_bright_l(ii);
    if exposure_time_l_tmp>1/30
        exposure_time_l_tmp=1/30;
    elseif exposure_time_l_tmp<0.0001
        exposure_time_l_tmp=0.0001;
    end
else exposure_time_l_tmp=exposure_time_l;
end
if log_road_bright_r(ii)>150 || log_road_bright_r(ii)<100
    exposure_time_r_tmp=exposure_time_r*133/log_road_bright_r(ii);
    if exposure_time_r_tmp>1/30
        exposure_time_r_tmp=1/30;
    elseif exposure_time_r_tmp<0.0001
        exposure_time_r_tmp=0.0001;
    end
else exposure_time_r_tmp=exposure_time_r;
end
if (exposure_time_l~=exposure_time_l_tmp) || (exposure_time_r~=exposure_time_r_tmp)
    exposure_time_l=exposure_time_l_tmp;
    exposure_time_r=exposure_time_r_tmp;
    src1.Exposure=exposure_time_l;
    src2.Exposure=exposure_time_r;
    trigger([vid1 vid2]);
    leftimg=getdata(vid1);
    rightimg=getdata(vid2);
    log_abnormal(ii)=ii;
    imwrite(leftimg,strcat(save_dir,'\left',int2str(ii),'.jpg'),'jpg');
    imwrite(rightimg,strcat(save_dir,'\right',int2str(ii),'.jpg'),'jpg');
    leftimg=double(max(leftimg,[],3));
    rightimg=double(max(rightimg,[],3));
else
    imwrite(leftimg_tmp,strcat(save_dir,'\left',int2str(ii),'.jpg'),'jpg');
    imwrite(rightimg_tmp,strcat(save_dir,'\right',int2str(ii),'.jpg'),'jpg');
    log_abnormal(ii)=0;
end
log_exposure_time_l(ii)=exposure_time_l;
log_exposure_time_r(ii)=exposure_time_r;