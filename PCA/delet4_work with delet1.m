load 'C:\Users\a0107257\Desktop\arrow_orti.mat';
% arrow=[zeros(size(arrow,1),20) arrow zeros(size(arrow,1),20)];
[r,c]=find(arrow==1);

% % cen_r_ori=round(size(arrow,1)/2);
% % cen_c_ori=round(size(arrow,2)/2);
cen_r_ori=round(mean(r));
cen_c_ori=round(mean(c));
cen_ind=find((r==cen_r_ori & c==cen_c_ori));
rot_cen_r=r(cen_ind);
rot_cen_c=c(cen_ind);
theta=(-45:1:45)/180*pi;
img=zeros(size(arrow,1),size(arrow,2),length(theta));
img_temp=zeros(63,82,length(theta));
img_pca=zeros(63*82,length(theta));
counter=0;
for i=1:length(theta)
    R=[cos(theta(i)) sin(theta(i)); -sin(theta(i)) cos(theta(i))];
    rot_rc=round(R*[r';c']);
    rot_r=rot_rc(1,:);
    rot_c=rot_rc(2,:);
    del_r=rot_r(cen_ind)-rot_cen_r;
    del_c=rot_c(cen_ind)-rot_cen_c;
    rot_r=rot_r-del_r;
    rot_c=rot_c-del_c;
    indtmp=(rot_r<=0 | rot_r>size(arrow,1) | rot_c<=0 | rot_c>size(arrow,2));
    rot_r(indtmp)=[];
    rot_c(indtmp)=[];
    indtmp=sub2ind(size(arrow),rot_r,rot_c);
    imgtmp=zeros(size(arrow));
    imgtmp(indtmp)=1;
    imgtmp=imfill(imgtmp,'holes');
    img(:,:,i)=imgtmp;
    img_temp(:,:,i)=img(round(1:11.27:710),round(1:6.09:499),i);
    a=img_temp(:,:,i);
    counter=counter+1;
    img_pca(:,counter)=a(:);
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
    imshow(reshape(img_pca(:,counter),63,82),[]);
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