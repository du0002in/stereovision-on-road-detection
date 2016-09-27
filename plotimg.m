% % % % % leftimg3d=zeros(ImgSize1/2,ImgSize2,3);
% % % % % leftimg3d(:,:,1)=leftimg(1:2:ImgSize1,:);
% % % % % leftimg3d(:,:,2)=leftimg(1:2:ImgSize1,:);
% % % % % leftimg3d(:,:,3)=leftimg(1:2:ImgSize1,:);
% % % % % leftimg3d=uint8(leftimg3d);
% % % % % xx=-CorepxL+cc_x_left;
% % % % % if leftrightlane==0
% % % % %     yy=(sel_ytmpr+1)/2;
% % % % % else
% % % % %     yy=(sel_ytmpl+1)/2;
% % % % % end
% % % % % ind=uint64(sub2ind([ImgSize1/2, ImgSize2],-yy,round(xx)));
% % % % % ind2=uint64(sub2ind([ImgSize1/2, ImgSize2],[-ytmpl -ytmpl2],round([xtmpl xtmpl2])));
% % % % % % ind3=sub2ind([ImgSize1/2, ImgSize2],-tmpvl,tmpul);
% % % % % leftimg3d(ind2)=255;
% % % % % leftimg3d(ind2+ImgSize1/2*ImgSize2)=0;
% % % % % leftimg3d(ind2+ImgSize1/2*ImgSize2*2)=0;
% % % % % % leftimg3d(ind3)=0;
% % % % % % leftimg3d(ind3+ImgSize1/2*ImgSize2)=255;
% % % % % % leftimg3d(ind3+ImgSize1/2*ImgSize2*2)=0;
% % % % % leftimg3d(ind)=0;
% % % % % leftimg3d(ind+ImgSize1/2*ImgSize2)=0;
% % % % % leftimg3d(ind+ImgSize1/2*ImgSize2*2)=0;
% % % % % % imwrite(leftimg3d,strcat('D:\Stereo\Image Collection\FresultL',int2str(ii),'.jpg'),'jpg');
% % % % % % h=subplot(2,1,1);imshow(leftimg3d,[]);
% % % % % 
% % % % % rightimg3d=zeros(ImgSize1/2,ImgSize2,3);
% % % % % rightimg3d(:,:,1)=rightimg(1:2:ImgSize1,:);
% % % % % rightimg3d(:,:,2)=rightimg(1:2:ImgSize1,:);
% % % % % rightimg3d(:,:,3)=rightimg(1:2:ImgSize1,:);
% % % % % rightimg3d=uint8(rightimg3d);
% % % % % xx=-CorepxR+cc_x_right;
% % % % % if leftrightlane==0
% % % % %     yy=(sel_ytmpr+1)/2;
% % % % % else
% % % % %     yy=(sel_ytmpl+1)/2;
% % % % % end
% % % % % ind=uint64(sub2ind([ImgSize1/2, ImgSize2],-yy,round(xx)));
% % % % % ind2=uint64(sub2ind([ImgSize1/2, ImgSize2],[-ytmpr -ytmpr2],round([xtmpr xtmpr2])));
% % % % % % ind3=sub2ind([ImgSize1/2, ImgSize2],-tmpvr,tmpur);
% % % % % rightimg3d(ind2)=255;
% % % % % rightimg3d(ind2+ImgSize1/2*ImgSize2)=0;
% % % % % rightimg3d(ind2+ImgSize1/2*ImgSize2*2)=0;
% % % % % % rightimg3d(ind3)=0;
% % % % % % rightimg3d(ind3+ImgSize1/2*ImgSize2)=255;
% % % % % % rightimg3d(ind3+ImgSize1/2*ImgSize2*2)=0;
% % % % % rightimg3d(ind)=0;
% % % % % rightimg3d(ind+ImgSize1/2*ImgSize2)=0;
% % % % % rightimg3d(ind+ImgSize1/2*ImgSize2*2)=0;
% % % % % % imwrite(leftimg3d,strcat('D:\Stereo\Image Collection\FresultR',int2str(ii),'.jpg'),'jpg');
% % % % % % subplot(2,1,2);imshow(rightimg3d,[]);
% % % % % % saveas(h,strcat(save_dir,'\Fresult',int2str(ii),'.jpg'),'jpg');
% % % % % % clf;
if miss_infor_flag==0
    xx=-CorepxR+cc_x_right;
    if leftrightlane==0
        yy=(sel_ytmpr+1)/2;
    else
        yy=(sel_ytmpl+1)/2;
    end
end

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
% % % rimgtmp3d(ind2)=255;
% % % rimgtmp3d(ind2+ImgSize1/2*ImgSize2)=0;
% % % rimgtmp3d(ind2+ImgSize1/2*ImgSize2*2)=0;
% % % rimgtmp3d(ind)=0;5
% % % rimgtmp3d(ind+ImgSize1/2*ImgSize2)=0;
% % % rimgtmp3d(ind+ImgSize1/2*ImgSize2*2)=0;
clf;
imshow(rimgtmp3d,[]);
if miss_infor_flag==0
    hold on;
    % % plot(inlier_xy_r(:,1),-inlier_xy_r(:,2),'k.');
    plot(xtmpr(ytmpr<-40),-ytmpr(ytmpr<-40),'g','LineWidth',2);
    plot(xtmpr2(ytmpr2<-40),-ytmpr2(ytmpr2<-40),'g','LineWidth',2);
    plot(xx, -yy, 'b.');
    if exist('indtmp_dumy','var')~=0
    indtmp2=find(indtmp_dumy==1);
    i_dumy=ir(indtmp2(find(Zw<7,1)));
    plot(xtmpr(i_dumy),-ytmpr(i_dumy),'r*');
    end
    hold off;
end
if ~isempty(ximg)
    if ~isempty(log_xyimg{ii})
        hold on;plot(ximg, yimg/2,'b','LineWidth',2);
        text(mean(ximg(1:4)), mean(yimg(1:4))/2,log_atrribute,'BackgroundColor',[.7 .9 .7]);
        hold off;
    else
        hold on;plot(ximg, yimg/2,'b--','LineWidth',2);
        text(mean(ximg(1:4)), mean([yimg(1:2) 480 480])/2,[log_data_name{end} ' (Estimated)'],'BackgroundColor',[.7 .9 .7]);
        hold off;
    end
    hold on;text(430,-11,'Warning Letter','BackgroundColor',[1 .0 .0]);hold off;
    if isempty(find(letter_gt==ii, 1))
        FP_letter=FP_letter+1;
    else
        TP_letter=TP_letter+1;
    end
else
    hold on;
    text(430,-11,'Warning Letter','BackgroundColor',[.8 .80 .80],'Color',[.7 .7 .7],'Edge',.7*[1 1 1]);
    hold off;
    if isempty(find(letter_gt==ii, 1))
        TN_letter=TN_letter+1;
    else
        FN_letter=FN_letter+1;
    end
end

if ~isempty(rr_arrow)
    hold on;plot(ximg_arrow,yimg_arrow/2,'r','LineWidth',2);
    text(605,-11,'Arrow','BackgroundColor',[1 .0 .0]);
    hold off;
    if isempty(find(arrow_gt==ii, 1))
        FP_arrow=FP_arrow+1;
    else
        TP_arrow=TP_arrow+1;
    end
    log_arrow(ii)=1;
elseif isempty(rr_arrow) && ~isempty(rr_tip)
    hold on;plot(ximg_arrow_tip,yimg_arrow_tip/2,'r','LineWidth',2);
    text(605,-11,'Arrow','BackgroundColor',[1 .0 .0]);
    hold off;
    if isempty(find(arrow_gt==ii, 1))
        FP_arrow=FP_arrow+1;
    else
        TP_arrow=TP_arrow+1;
    end
    log_arrow(ii)=1;
elseif isempty(rr_arrow) && isempty(rr_tip) && log_p_tip(ii-1)>0
    hold on;plot(ximg_arrow_tip_est,yimg_arrow_tip_est/2,'r--','LineWidth',2);
    text(605,-11,'Arrow','BackgroundColor',[1 .0 .0]);
    hold off;
    if isempty(find(arrow_gt==ii, 1))
        FP_arrow=FP_arrow+1;
    else
        TP_arrow=TP_arrow+1;
    end
    log_arrow(ii)=1;
else
    hold on;
    text(605,-11,'Arrow','BackgroundColor',[.8 .80 .80],'Color',[.7 .7 .7],'Edge',.7*[1 1 1]);
    hold off;
    if isempty(find(arrow_gt==ii, 1))
        TN_arrow=TN_arrow+1;
    else
        FN_arrow=FN_arrow+1;
    end
    log_arrow(ii)=0;
end

if log_zebra(ii)==1
    hold on;plot(ximg_zebra, yimg_zebra/2, 'y', 'LineWidth', 2);
    text(566,-11,'Zebra','BackgroundColor',[1 .0 .0]);
    hold off;
    if isempty(find(zebra_gt==ii, 1))
        FP_zebra=FP_zebra+1;
    else
        TP_zebra=TP_zebra+1;
    end
else
    hold on;
    text(566,-11,'Zebra','BackgroundColor',[.8 .80 .80],'Color',[.7 .7 .7],'Edge',.7*[1 1 1]);
    hold off;
    if isempty(find(zebra_gt==ii, 1))
        TN_zebra=TN_zebra+1;
    else
        FN_zebra=FN_zebra+1;
    end
end

if log_hump(ii)==1
    hold on;plot(ximg_hump, yimg_hump/2, 'c', 'LineWidth', 2); 
    text(524,-11,'Hump','BackgroundColor',[1 .0 .0]);
    hold off;
    if isempty(find(hump_gt==ii, 1))
        FP_hump=FP_hump+1;
    else
        TP_hump=TP_hump+1;
    end
else
    hold on;
    text(524,-11,'Hump','BackgroundColor',[.8 .80 .80],'Color',[.7 .7 .7],'Edge',.7*[1 1 1]);
    hold off;
    if isempty(find(hump_gt==ii, 1))
        TN_hump=TN_hump+1;
    else
        FN_hump=FN_hump+1;
    end
end

hold on;text(1, -11, ['Dist to right: ' num2str(log_Estimate_xc(ii),'%1.2f') 'm'],'BackgroundColor',[.7 .9 .7]);
text(125, -11, ['Moving direction: ' num2str(log_Estimate_d_theta(ii),'%1.3f') 'rad'],'BackgroundColor',[.7 .9 .7]);
hold off;

% % % % if length(hump_row)>=20
% % % % %     hold on;plot([1 ImgSize2 ImgSize2 1 1],[hump_row(1) hump_row(1) hump_row(end) hump_row(end) hump_row(1)],'r');
% % % %     hold on;
% % % %     for hump_i=[1:length(hump_row)]
% % % %         plot([1 640],[hump_row(hump_i) hump_row(hump_i)],'r');
% % % %     end
% % % %     hold off;
% % % % end
% % % hold on;plot(expos_box{3},-expos_box{7},'--');
% % % plot(expos_box{4},-expos_box{8},'--');hold off;
% % % hold on;plot([1 640],[25 25],'r--');
% % % plot([1 640],[55 55],'r--');
% % % plot([1 640],[120 120],'r--');hold off;
frame=getframe(gcf);
im=frame2im(frame);
imwrite(im,strcat(save_dir,'\3DFresult_letter',int2str(ii),'.jpg'),'jpg');
