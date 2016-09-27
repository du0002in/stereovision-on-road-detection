clear all;
load 'D:\StereoLaptop_new\PCA\arrow_orti.mat';
figure;
% arrow=[zeros(size(arrow,1),20) arrow zeros(size(arrow,1),20)];
arrow(:,1:72)=[]; arrow(:,(end-71):end)=[];
arrow(296:end,:)=[];
arrow(:,[1:91 (end-91):end])=[];
vetex=[1 710 710 1 1;52 52 120 120 52];
% % vetex=[1 295 295 1 1;52 52 120 120 52]; %for arrow tip only

[r,c]=find(arrow==1);

% % cen_r_ori=round(size(arrow,1)/2);
% % cen_c_ori=round(size(arrow,2)/2);
cen_r_ori=round(mean(r));
cen_c_ori=round(mean(c));
cen_ind=find((r==cen_r_ori & c==cen_c_ori));
rot_cen_r=r(cen_ind);
rot_cen_c=c(cen_ind);
theta=(-30:2:30)/180*pi;
img=zeros(size(arrow,1),size(arrow,2),length(theta));

row_res=round(1:(1*11.27):295); 
col_res=round(1:(1*6.09):172);
img_temp=zeros(length(row_res),length(col_res),length(theta));
img_pca=zeros(length(row_res)*length(col_res),length(theta));
% % img_temp=zeros(63,59,length(theta));
% % img_pca=zeros(63*59,length(theta));

counter=0;
for i=1:length(theta)
    R=[cos(theta(i)) sin(theta(i)); -sin(theta(i)) cos(theta(i))];
    rot_rc=round(R*[r';c']);
    vetex_rot=R*vetex;
    rot_r=rot_rc(1,:);
    rot_c=rot_rc(2,:);
    del_r=rot_r(cen_ind)-rot_cen_r;
    del_c=rot_c(cen_ind)-rot_cen_c;
    rot_r=rot_r-del_r;
    rot_c=rot_c-del_c;
    vetex_rot(1,:)=vetex_rot(1,:)-del_r;
    vetex_rot(2,:)=vetex_rot(2,:)-del_c;
    indtmp=(rot_r<=0 | rot_r>size(arrow,1) | rot_c<=0 | rot_c>size(arrow,2));
    rot_r(indtmp)=[];
    rot_c(indtmp)=[];
    indtmp=sub2ind(size(arrow),rot_r,rot_c);
    imgtmp=zeros(size(arrow));
    imgtmp(indtmp)=1;
    imgtmp=imfill(imgtmp,'holes');
    img(:,:,i)=imgtmp;
    img_temp(:,:,i)=img(row_res,col_res,i);
    vetex_rot(1,:)=(vetex_rot(1,:)-1)/1.0/11.27+1;
    vetex_rot(2,:)=(vetex_rot(2,:)-1)/1.0/6.09+1;
    a=img_temp(:,:,i);
    counter=counter+1;
    img_pca(:,counter)=a(:);
    log_a(counter)=sum(a(:));
    sum_a_col=sum(a);
% % % %     for j=1:length(sum_a_col)
% % % %         if sum_a_col(j)==0
% % % %             atmp=[a(:,j+1:end) zeros(size(a,1),j)];
% % % %             counter=counter+1;
% % % %             img_pca(:,counter)=atmp(:);
% % % %         else
% % % %             break;
% % % %         end
% % % %     end
% % % %     for j=length(sum_a_col):-1:1
% % % %         if sum_a_col(j)==0
% % % %             atmp=[zeros(size(a,1),size(a,2)-j+1), a(:,1:(j-1))];
% % % %             counter=counter+1;
% % % %             img_pca(:,counter)=atmp(:);
% % % %         else
% % % %             break;
% % % %         end
% % % %     end
    imshow(reshape(img_pca(:,counter),length(row_res),length(col_res)),[]);
    hold on;plot(vetex_rot(2,:), vetex_rot(1,:),'r');hold off;
    log_vetex{counter}=vetex_rot;
    drawnow; pause(0.1);
end
img_pca_mean=mean(img_pca,2);
X=zeros(size(img_pca));
for i=1:length(theta)
    X(:,i)=img_pca(:,i)-img_pca_mean;
end
[u,s,v]=svd(X,'econ');
g=s*v';
% [coe,sco,lat]=pca(X);