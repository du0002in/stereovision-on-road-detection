% % maxconvimg=dumy1;
% % if ~isempty(xtmpl2) && ~isempty(xtmpr2)
% %     xtmpl2=[xtmpl2 ImgSize2*ones(1,ImgSize1/2+ytmpl2(end))];
% %     xtmpr2=[xtmpr2 ImgSize2*ones(1,ImgSize1/2+ytmpr2(end))];
% %     maxconvimg=C_Search_region(maxconvimg,xtmpl2,xtmpr2+ImgSize2+30,-ytmpr2(1),ImgSize1/2);
% % else
% %     if leftrightlane==0
% %         xtmpl_dumy=expos_box{1}; xtmpr_dumy=expos_box{3};
% %         maxconvimg=C_Search_region(maxconvimg,xtmpl_dumy,xtmpr_dumy+ImgSize2+30,-expos_box{8}(1),ImgSize1/2);
% %     else
% %         xtmpl_dumy=expos_box{2}; xtmpr_dumy=expos_box{4};
% %         maxconvimg=C_Search_region(maxconvimg,xtmpl_dumy,xtmpr_dumy+ImgSize2+30,-expos_box{7}(1),ImgSize1/2);
% %     end
% % end
% % num_ite=1960;
% % [sel_xtmpl,sel_ytmpl, xtmpl, ytmpl,sel_xtmpr,sel_ytmpr, xtmpr, ytmpr, num_in, para, inlier_xy_r, paraz,break_location,log_dumy(ii,:)]= ...
% %     linefitting_simu2(maxconvimg(:,1:ImgSize2),maxconvimg(:,end-ImgSize2+1:end), D_pre, H, num_ite,pre_line_type,logC0(ii-1));
% % if num_in>=60
% %     return;
% % end

if mod(floor(D_pre),2)==0
    D_pre1=floor(D_pre)-1;
else D_pre1=floor(D_pre);
end
yy=min(floor(D_pre1),top_row):-1:btm_row;
if marking_width==0.1
    marking_width1=0.15;
else marking_width1=0.1;
end
tdxtmp1=(cos(Phai)/H*(D_pre-yy)*marking_width1)/2; %0.10m is the lane marking width
if length(tdxtmp1)>=length(tdx1)
    tdx1=round(tdxtmp1((end-length(tdx1)+1):end));
else
    tdx1=round([0.5*ones(1,(length(tdx1)-length(tdxtmp1))) tdxtmp1]);
end
tdx1(tdx1==0)=1;
[~,ix1,~]=unique(tdx1,'stable');ix1=ix1';
j=1;
for i=ix1
    weight3=exp((-0.5*(-tdx1(i):tdx1(i)).^2/(tdx1(i))^2));
    weight3=weight3/sum(weight3);
    xfil3{j}=weight3;
    j=j+1;
end
ix1(length(ix1)+1)=ImgSize1-Offset+1;
ix1=round(ix1/2)+Offset/2;

ImgSmooth=zeros(size(img));
for i=1:(length(ix1)-1)
    ImgSmooth(ix1(i):(ix1(i+1)-1),:)=conv2(img(ix1(i):(ix1(i+1)-1),:),xfil3{i},'same');
end
maxconvimg1=C_new_Ridge(ImgSmooth,tdx1,Offset);
maxconvimg1=maxconvimg1.*mask;
C_refine_Ridge(maxconvimg1,Offset);
if ~isempty(ximg)
    ximg1=ximg(1:4)+ImgSize2+30; yimg1=yimg(1:4)/2;
    min_x=max([ImgSize2+31, min(ximg1)]);   max_x=min([2*ImgSize2+30, max(ximg1)]);
    min_y=max([1, min(yimg1)]);             max_y=min([ImgSize1/2, max(yimg1)]);
    maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
    maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
    [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
    if ~isempty(maxconvimgtmp_x)
        maxconvimgtmp_y=maxconvimgtmp_y-1;
        maxconvimgtmp_x=maxconvimgtmp_x-1;
        ximg1=ximg1-1;            yimg1=yimg1-1;
        C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
    end
    log_letter(ii)=1;
end

dumy2=maxconvimg1(1:end,1:end); %after letter detection

% log_zebra(ii)=0;
if ~(log_hump(ii-1)==1 || log_hump(ii-2)==1)
%     zebra_recognition;
    if log_zebra(ii)==1
        ximg1=ximg_zebra(1:4)+ImgSize2+30; yimg1=yimg_zebra(1:4)/2;
        min_x=max([ImgSize2+31, min(ximg1)]);   max_x=min([2*ImgSize2+30, max(ximg1)]);
        min_y=max([1, min(yimg1)]);             max_y=min([ImgSize1/2, max(yimg1)]);
        maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
        maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
        [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
        if ~isempty(maxconvimgtmp_x)
            maxconvimgtmp_y=maxconvimgtmp_y-1;
            maxconvimgtmp_x=maxconvimgtmp_x-1;
            ximg1=ximg1-1;            yimg1=yimg1-1;
            C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
        end
    end
end

% log_hump(ii)=0;
if ~(log_zebra(ii-1)==1 || log_zebra(ii)==1)
%     hump_recognition;
    if log_hump(ii)==1
        ximg1=ximg_hump(1:4)+ImgSize2+30; yimg1=yimg_hump(1:4)/2;
        min_x=max([ImgSize2+31, min(ximg1)]);   max_x=min([2*ImgSize2+30, max(ximg1)]);
        min_y=max([1, min(yimg1)]);             max_y=min([ImgSize1/2, max(yimg1)]);
        maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
        maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
        [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
        if ~isempty(maxconvimgtmp_x)
            maxconvimgtmp_y=maxconvimgtmp_y-1;
            maxconvimgtmp_x=maxconvimgtmp_x-1;
            ximg1=ximg1-1;            yimg1=yimg1-1;
            C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
        end
        ximg1=ximg_hump(1:4); yimg1=yimg_hump(1:4)/2;
        min_x=max([1, min(ximg1)]);   max_x=min([ImgSize2, max(ximg1)]);
        min_y=max([1, min(yimg1)]);   max_y=min([ImgSize1/2, max(yimg1)]);
        maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
        maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
        [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
        if ~isempty(maxconvimgtmp_x)
            maxconvimgtmp_y=maxconvimgtmp_y-1;
            maxconvimgtmp_x=maxconvimgtmp_x-1;
            ximg1=ximg1-1;            yimg1=yimg1-1;
            C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
        end
    end
end

% arrow_recognition;
if ~isempty(rr_arrow)
    ximg1=ximg_arrow(1:4)+ImgSize2+30; yimg1=yimg_arrow(1:4)/2;
    min_x=max([ImgSize2+31, min(ximg1)]);   max_x=min([2*ImgSize2+30, max(ximg1)]);
    min_y=max([1, min(yimg1)]);             max_y=min([ImgSize1/2, max(yimg1)]);
    maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
    maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
    [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
    if ~isempty(maxconvimgtmp_x)
        maxconvimgtmp_y=maxconvimgtmp_y-1;
        maxconvimgtmp_x=maxconvimgtmp_x-1;
        ximg1=ximg1-1;            yimg1=yimg1-1;
        C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
    end
elseif isempty(rr_arrow) && ~isempty(rr_tip)
    ximg1=ximg_arrow_tip(1:4)+ImgSize2+30; yimg1=yimg_arrow_tip(1:4)/2;
    min_x=max([ImgSize2+31, min(ximg1)]);   max_x=min([2*ImgSize2+30, max(ximg1)]);
    min_y=max([1, min(yimg1)]);             max_y=min([ImgSize1/2, max(yimg1)]);
    maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
    maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
    [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
    if ~isempty(maxconvimgtmp_x)
        maxconvimgtmp_y=maxconvimgtmp_y-1;
        maxconvimgtmp_x=maxconvimgtmp_x-1;
        ximg1=ximg1-1;            yimg1=yimg1-1;
        C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
    end
elseif isempty(rr_arrow) && isempty(rr_tip) && log_p_tip(ii-1)>0
    ximg1=ximg_arrow_tip_est(1:4)+ImgSize2+30; yimg1=yimg_arrow_tip_est(1:4)/2;
    min_x=max([ImgSize2+31, min(ximg1)]);   max_x=min([2*ImgSize2+30, max(ximg1)]);
    min_y=max([1, min(yimg1)]);             max_y=min([ImgSize1/2, max(yimg1)]);
    maxconvimgtmp=zeros(ImgSize1/2,2*ImgSize2+30);
    maxconvimgtmp(min_y:max_y,min_x:max_x)=maxconvimg1(min_y:max_y,min_x:max_x);
    [maxconvimgtmp_y, maxconvimgtmp_x]=find(maxconvimgtmp==1);
    if ~isempty(maxconvimgtmp_x)
        maxconvimgtmp_y=maxconvimgtmp_y-1;
        maxconvimgtmp_x=maxconvimgtmp_x-1;
        ximg1=ximg1-1;            yimg1=yimg1-1;
        C_in_poly_new(maxconvimg1, maxconvimgtmp_x, maxconvimgtmp_y, ximg1, yimg1);
    end
end
dumy5=maxconvimg1(1:end,1:end);
dumy1=maxconvimg1;

num_ite=15000; %number of iterations for RANSAC
if ~isempty(xtmpl) && ~isempty(xtmpr)
    xtmpl=[xtmpl ImgSize2*ones(1,ImgSize1/2+ytmpl(end))];
    xtmpr=[xtmpr ImgSize2*ones(1,ImgSize1/2+ytmpr(end))];
    maxconvimg1=C_Search_region(maxconvimg1,xtmpl,xtmpr+ImgSize2+30,-ytmpr(1),ImgSize1/2);
    %         maxconvimg=C_Search_region(maxconvimg,xtmpl,xtmpr+ImgSize2+30,-ytmpr(1),-ytmpr(end));
    %         maxconvimg(-ytmpr(end):end,:)=0;
    maxconvimg1(1:-ytmpr(1),:)=0;
    maxconvimg1(-ytmpr(end):ImgSize1/2,:)=0;
    num_ite=960;
else
    if leftrightlane==0
        xtmpl_dumy=expos_box{1}; xtmpr_dumy=expos_box{3};
        maxconvimg1=C_Search_region(maxconvimg1,xtmpl_dumy,xtmpr_dumy+ImgSize2+30,-expos_box{7}(1),ImgSize1/2);
    else
        xtmpl_dumy=expos_box{2}; xtmpr_dumy=expos_box{4};
        maxconvimg1=C_Search_region(maxconvimg1,xtmpl_dumy,xtmpr_dumy+ImgSize2+30,-expos_box{8}(1),ImgSize1/2);
    end
end
dumy4=maxconvimg1; %after removing lane center pixels based on pre reseults, then combine with pre dumys
dumy=[dumy1;dumy2];
imwrite([dumy3;dumy2;dumy5;dumy1;dumy4],strcat(save_dir,'\oriRidge',int2str(ii),'.jpg'),'jpg');
[sel_xtmpl,sel_ytmpl, xtmpl, ytmpl,sel_xtmpr,sel_ytmpr, xtmpr, ytmpr, num_in, para, inlier_xy_r, paraz,break_location, ~]= ...
    linefitting_simu2(maxconvimg1(:,1:ImgSize2),maxconvimg1(:,end-ImgSize2+1:end), D_pre, H, num_ite,pre_line_type,logC0(ii-1),log_zebra(ii));
% lognum_in(ii)=num_in;
if num_in>=60
    tdx=tdx1;
    xfil2=xfil3;
    ix=ix1;
    marking_width=marking_width1;
    return;
end


 maxconvimg1=dumy1;

if ~isempty(xtmpl2) && ~isempty(xtmpr2)
    xtmpl2=[xtmpl2 ImgSize2*ones(1,ImgSize1/2+ytmpl2(end))];
    xtmpr2=[xtmpr2 ImgSize2*ones(1,ImgSize1/2+ytmpr2(end))];
    maxconvimg1=C_Search_region(maxconvimg1,xtmpl2,xtmpr2+ImgSize2+30,-ytmpr2(1),ImgSize1/2);
else
    if leftrightlane==0
        xtmpl_dumy=expos_box{2}; xtmpr_dumy=expos_box{4};
        maxconvimg1=C_Search_region(maxconvimg1,xtmpl_dumy,xtmpr_dumy+ImgSize2+30,-expos_box{8}(1),ImgSize1/2);
    else
        xtmpl_dumy=expos_box{1}; xtmpr_dumy=expos_box{3};
        maxconvimg1=C_Search_region(maxconvimg1,xtmpl_dumy,xtmpr_dumy+ImgSize2+30,-expos_box{7}(1),ImgSize1/2);
    end
end
num_ite=1960;
[sel_xtmpl,sel_ytmpl, xtmpl, ytmpl,sel_xtmpr,sel_ytmpr, xtmpr, ytmpr, num_in, para, inlier_xy_r, paraz,break_location,~]= ...
            linefitting_simu2(maxconvimg1(:,1:ImgSize2),maxconvimg1(:,end-ImgSize2+1:end), D_pre, H, num_ite,pre_line_type,logC0(ii-1),log_zebra(ii));

if num_in>=60
    if leftrightlane==0
        leftrightlane=1;
    else
        leftrightlane=0;
    end
    tdx=tdx1;
    xfil2=xfil3;
    ix=ix1;
    marking_width=marking_width1;
end
