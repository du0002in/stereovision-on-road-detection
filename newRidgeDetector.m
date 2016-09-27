clear all;
warning off;

% leftimg_dir='D:\Stereo\Image Collection4\left';
% rightimg_dir='D:\Stereo\Image Collection4\right';
% save_dir='D:\Stereo\Image Collection4';
leftimg_dir='D:\Stereo\Image Collection6\right';
rightimg_dir='D:\Stereo\Image Collection6\left';
save_dir='D:\Stereo\Image Collection6\test';
% leftimg_dir='D:\Stereo\realrun\test4\left';
% rightimg_dir='D:\Stereo\realrun\test4\right';
% save_dir='D:\Stereo\realrun\test4\postrun';
% leftimg_dir='D:\Stereo\Calibration\realcali\realleft';
% rightimg_dir='D:\Stereo\Calibration\realcali\realright';
% save_dir='D:\Stereo\Calibration\realcali';
figure;
%for collection6 ground truth
log_ite_s=[575 695 1110 1490];
log_ite_e=[665 960 1444 1540];
Offset=40;
% log_ite_s=1;log_ite_e=5;
%for collection6
% log_ite_s=[346    555 1110 1725];
% log_ite_e=[382  960 1545 2040];

weight2y=exp((-0.5*(-1:1)'.^2/1^2));
weight2y=weight2y/sum(weight2y);
yfil2=weight2y;

for outloop=1:1
    ite_s=log_ite_s(outloop);
    ite_e=log_ite_e(outloop);
    global_varibles;
    InitializationGlobalVariable;
    load 'cali3.mat';
    % load(strcat(save_dir,'\odemotry.mat'));%actural data readings from EV velocity and steering
    load(strcat(save_dir,'\log_tdx.mat'));
    
    D_pre=852*tan(.2757);H=1.665;marking_width_pre=0.15;diff_marking_width=0;
    marking_width_pre2=0.15;continue_frame=0;
    top_row=2*floor((cc_y-40)/2)-1; %make sure cc_y-40 is odd
    btm_row=2*ceil((cc_y-480)/2); %make sure cc_y-480 is even
    cc_y=2*ceil(cc_y/2);
    % D_pre=174.4855; H=1.6428;
    SampleRectification;
    Ridge_width_ini=8.5;
    ImgSize1=480; ImgSize2=640;
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
    
    for ii=ite_s:ite_e
        if ~isempty(log_tdx{ii})
            tdx=log_tdx{ii};
        end
        [~,ix,~]=unique(tdx,'stable');ix=ix';
        j=1;
        for i=ix
            weight3=exp((-0.5*(-tdx(i):tdx(i)).^2/(tdx(i))^2));
            weight3=weight3/sum(weight3);
            xfil2{j}=weight3;
            j=j+1;
        end
        ix(length(ix)+1)=ImgSize1-Offset+1;
        ix=round(ix/2)+Offset/2;
        leftimg=imread(strcat(leftimg_dir,int2str(ii),'.jpg'));
        rightimg=imread(strcat(rightimg_dir,int2str(ii),'.jpg'));
        leftimg=double(max(leftimg,[],3));
        rightimg=double(max(rightimg,[],3));
        % tic
        leftimg=Retify(leftimg,a1_left,a2_left,a3_left,a4_left,ind_1_left,ind_2_left,ind_3_left,ind_4_left,ind_new_left);
        rightimg=Retify(rightimg,a1_right,a2_right,a3_right,a4_right,ind_1_right,ind_2_right,ind_3_right,ind_4_right,ind_new_right);
        img=[leftimg(1:2:ImgSize1,:) zeros(ImgSize1/2, 30) rightimg(1:2:ImgSize1,:)];
        img(235:end,:)=0;
        ImgSmooth=zeros(size(img));
        for i=1:(length(ix)-1)
            ImgSmooth(ix(i):(ix(i+1)-1),:)=conv2(img(ix(i):(ix(i+1)-1),:),xfil2{i},'same');
        end
        tic;
        Fimg=C_new_Ridge(ImgSmooth,tdx,Offset);
        logt(ii)=toc;
        tic;
        ImgSmooth=conv2(ImgSmooth,yfil2,'same');
        [kcsditmp,th]=C_Ridge(ImgSmooth,306400);
        logt2(ii)=toc;
%         Fimg=zeros(size(ImgSmooth));
%         for i=(Offset/2+1):234
%             if i==92
%                 a=0;
%             end
%             aa=ImgSmooth(i,:);
%             width=2*tdx(2*i-Offset);
%             for j=(width+1):(length(aa)-width-1)
%                 if j==133
%                     a=0;
%                 end
%                 rg=aa(j-width:j+width);
%                 [~,ind]=max(rg);
%                 if ind~=(width+1)
%                     continue;
%                 else
%                     ll=rg(2:ind)-rg(1:ind-1);
%                     rr=rg(ind:end-1)-rg(ind+1:end);
%                     ll((ll<3) & (ll>-3))=0;
%                     rr((rr<3) & (rr>-3))=0;
%                     if (sum(ll<0)<=0) && (sum(rr<0)<=0) && (rg(ind)-min(rg)>30) && (abs(rg(1)-rg(end))<=25)
%                         Fimg(i,j)=1;
%                     end
%                 end
%             end
%         end
        imwrite(Fimg,strcat(save_dir,'\Fimg_new',int2str(ii),'.jpg'),'jpg');
        imshow(Fimg,[]);
        pause(0.002);
    end
end
