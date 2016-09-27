clear all;
% matlabpool('open',4);
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
load cali3.mat;
% load 'all_char_at_all_angle.mat';
load 'roadmarkingtemplates_low.mat';
SampleRectification;
f='D:\StereoLaptop_new\Image Collection(new Cam)8';
save_dir=f;
figure;
% f='D:\Stereo (Laptop)\RealRun_Laptop_new_Cam\RealRunImg\test166';
rightimg_dir=[f '\right'];
[ximg,yimg]=meshgrid(1:640,1:2:480);
xxc=-ximg+cc_x_right;
yyc=-yimg+cc_y;
total_counter=0;
ximg=[]; yimg=[];
% for ii=2:362
for ii=460:772
% for ii=118:526
% for ii=702:1085
%     for ii=348:350

    rightimg=imread(strcat(rightimg_dir,int2str(ii),'.jpg'));
    % figure;imshow(rightimg,[]);
    rightimg=double(max(rightimg,[],3));
    rightimg=Retify(rightimg,a1_right,a2_right,a3_right,a4_right,ind_1_right,ind_2_right,ind_3_right,ind_4_right,ind_new_right);
    % ii=ii-1;
    % p=logPhai_f(ii); H=logH_f(ii); focal=684; d_theta=log_Estimate_d_theta(ii);
    temp_p=0.3; temp_H=1.66; temp_td=230;
    template_matching;
    imshow(tempimg,[]);
    if ~isempty(ximg)
        hold on;plot(ximg, yimg/2,'r');
        text(mean(ximg(1:4)), mean(yimg(1:4))/2,log_atrribute,'BackgroundColor',[.7 .9 .7]);
        hold off;
        pause(0.1);
    end
    if (isempty(ximg) && ~isempty(log_xyimg{ii-1}))
        if max(log_xyimg{ii-1}(2,:))>300
            letter_y_down_l=(log_xyimg{ii-1}(2,1)+log_xyimg{ii-1}(2,4))*0.5;
            letter_y_down_r=(log_xyimg{ii-1}(2,2)+log_xyimg{ii-1}(2,3))*0.5;
            letter_x_l=(log_xyimg{ii-1}(1,1)+log_xyimg{ii-1}(1,4))*0.5;
            letter_x_r=(log_xyimg{ii-1}(1,2)+log_xyimg{ii-1}(1,3))*0.5;
            esti_ximg(1)=letter_x_l;
            esti_ximg(2)=letter_x_r;
            esti_yimg(1)=letter_y_down_l;
            esti_yimg(2)=letter_y_down_r;
            esti_yimg(3)=max([log_xyimg{ii-1}(2,3)+letter_y_down_l-log_xyimg{ii-1}(2,1) 480]);
            esti_yimg(4)=max([log_xyimg{ii-1}(2,4)+letter_y_down_r-log_xyimg{ii-1}(2,2) 480]);
            esti_yimg(5)=esti_yimg(1);
            esti_ximg(3)=(esti_yimg(3)-log_xyimg{ii-1}(2,3))/(log_xyimg{ii-1}(2,3)-log_xyimg{ii-1}(2,2)) ...
                *(log_xyimg{ii-1}(1,3)-log_xyimg{ii-1}(1,2))+log_xyimg{ii-1}(1,3);
            esti_ximg(4)=(esti_yimg(4)-log_xyimg{ii-1}(2,4))/(log_xyimg{ii-1}(2,4)-log_xyimg{ii-1}(2,1)) ...
                *(log_xyimg{ii-1}(1,4)-log_xyimg{ii-1}(1,1))+log_xyimg{ii-1}(1,4);
            esti_ximg(5)=esti_ximg(1);
            ximg=esti_ximg; yimg=esti_yimg;
            hold on;plot(ximg, yimg/2, 'g');
            hold off; pause(0.1);
        end
    end
    frame=getframe(gcf);
    im=frame2im(frame);
    imwrite(im,strcat(save_dir,'\LetterMatch',int2str(ii),'.jpg'),'jpg');
    drawnow; 
end
% figure;plot(log_min); hold on;plot(log_min,'r*');
% set(gca,'XTick',[1:counter]);
% set(gca,'XTickLabel',log_atrribute);

% hold on;plot(ximg, yimg/2,'r');
% text(mean(ximg(1:4)), mean(yimg(1:4))/2,log_atrribute{counter},'BackgroundColor',[.7 .9 .7]);



