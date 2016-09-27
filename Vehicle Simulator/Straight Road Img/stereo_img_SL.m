%Generate the simulated stereo images given the vehicle global cor. and its
%direction w.r.t global cor. +x axis
function [img_left, img_right]=stereo_img_SL(X,Y,theta,phi,H,dist_left_right,L_wheel_cam,Tpix2cam_left, Tpix2cam_right, map_x,map_y,map_z, Lane_marking_point_x, Lane_marking_point_y, Lane_marking_point_z)

[Tcar2g,IDX,~]=Car2Gloable_SL(X,Y,theta, map_x, map_y, map_z);
[Tcam2car_left,Tcam2car_right]=Cam2Car(phi,H,dist_left_right,L_wheel_cam);
Ttmp_left=Tcam2car_left*Tcar2g;
Ttmp_right=Tcam2car_right*Tcar2g;
Tpix2g_left=Tpix2cam_left*Ttmp_left(1:3,:); %3x4
Tpix2g_right=Tpix2cam_right*Ttmp_right(1:3,:); %3x4

[pixX, pixY]=meshgrid(1:640,1:480);
pixX=pixX(:);pixY=pixY(:);

% load 'map_cordinate4_dotted';
FOV_deapth=3000; %the deapth of field of view. 600 data points~=30m
if (IDX+FOV_deapth)<=length(map_x)
    FOV_ind=IDX:(IDX+FOV_deapth);
else
    FOV_ind=[IDX:length(map_x) 1:(FOV_deapth-(length(map_x)-IDX))];
end
img_left=zeros(480,640);
img_right=zeros(480,640);
Xg=Lane_marking_point_x(FOV_ind,:);
Yg=Lane_marking_point_y(FOV_ind,:);
Zg=Lane_marking_point_z(FOV_ind,:);
Xpix_left=zeros(size(Xg,1),size(Xg,2));
Ypix_left=zeros(size(Xg,1),size(Xg,2));
% Zpix_left=zeros(size(Xg,1),size(Xg,2));
Xpix_right=zeros(size(Xg,1),size(Xg,2));
Ypix_right=zeros(size(Xg,1),size(Xg,2));


ii=[1:size(Xg,2)/2,size(Xg,2):-1:(size(Xg,2)/2+1)];
for itmp=1:size(Xg,2)
    i=ii(itmp);
    XYZpix_left=Tpix2g_left*[Xg(:,i) Yg(:,i) Zg(:,i) ones(size(Xg,1),1)]';
    Xpix_left(:,i)=round(XYZpix_left(1,:)'./XYZpix_left(3,:)');
    Ypix_left(:,i)=round(XYZpix_left(2,:)'./XYZpix_left(3,:)');
    indtmp1=(Xpix_left(:,i)>=1 & Xpix_left(:,i)<=640) & (Ypix_left(:,i)>=1 & Ypix_left(:,i)<=480);
    if mod(i,2)==0
        x_even=Xpix_left(:,i); y_even=Ypix_left(:,i);
        x_even(~indtmp1)=nan;
    else
        x_odd=Xpix_left(:,i); y_odd=Ypix_left(:,i);
        x_odd(~indtmp1)=nan;
    end
    if mod(i,2)==0
        array=~isnan(x_odd);
        edgeArray=diff([0; array; 0]);
        indices_odd=[find(edgeArray>0) find(edgeArray<0)-1];
        
        array=~isnan(x_even);
        edgeArray=diff([0;array;0]);
        indices_even=[find(edgeArray>0) find(edgeArray<0)-1];
        if ~isempty(indices_odd) && ~isempty(indices_even)
            for j=1:size(indices_odd,1)
                aa=indices_odd(j,1):indices_odd(j,2);
                for k=1:size(indices_even,1)
                    bb=indices_even(k,1):indices_even(k,2);
                    if ~((aa(end)<bb(1)) || (bb(end)<aa(1)))
                        
                        XV=[x_odd(aa);flipud(x_even(bb));x_odd(aa(1))];
                        YV=[y_odd(aa); flipud(y_even(bb));y_odd(aa(1))];
                        XV_max=max(XV); XV_min=min(XV);
                        YV_max=max(YV); YV_min=min(YV);
% %                         tic;
% %                         ind_x=((pixX<=XV_max) & (pixX>=XV_min));
% %                         ind_y=((pixY<=YV_max) & (pixY>=YV_min));
% %                         ind_xy=(ind_x & ind_y);
% %                         ptXtmp=pixX(ind_xy); ptYtmp=pixY(ind_xy);
% %                         toc
                        [ptXtmp,ptYtmp]=meshgrid(XV_min:XV_max,YV_min:YV_max);
                        ptXtmp=ptXtmp(:);ptYtmp=ptYtmp(:);
                        if issorted(ptYtmp)
                            [in_on,~]=inpoly([ptXtmp ptYtmp],[XV YV]);
                        elseif issorted(ptXtmp)
                            [in_on,~]=inpoly([ptYtmp ptXtmp],[YV XV]);
                        else
                            [in_on,~]=inpoly([ptXtmp ptYtmp],[XV YV]);
                        end
                        pix_in_on=sub2ind([480,640],[ptYtmp(in_on==1);YV],[ptXtmp(in_on==1);XV]);
                        img_left(pix_in_on)=1;

%                         patch([x_odd(aa);flipud(x_even(bb))], ...
%                             [y_odd(aa); flipud(y_even(bb))], 'w');
                        break;
                    end
                end
            end
        end
    end  
end

ii=[1:size(Xg,2)/2,size(Xg,2):-1:(size(Xg,2)/2+1)];
for itmp=1:size(Xg,2)
    i=ii(itmp);
    XYZpix_right=Tpix2g_right*[Xg(:,i) Yg(:,i) Zg(:,i) ones(size(Xg,1),1)]';
    Xpix_right(:,i)=round(XYZpix_right(1,:)'./XYZpix_right(3,:)');
    Ypix_right(:,i)=round(XYZpix_right(2,:)'./XYZpix_right(3,:)');
    indtmp1=(Xpix_right(:,i)>=1 & Xpix_right(:,i)<=640) & (Ypix_right(:,i)>=1 & Ypix_right(:,i)<=480);
    if mod(i,2)==0
        x_even=Xpix_right(:,i); y_even=Ypix_right(:,i);
        x_even(~indtmp1)=nan;
    else
        x_odd=Xpix_right(:,i); y_odd=Ypix_right(:,i);
        x_odd(~indtmp1)=nan;
    end
    if mod(i,2)==0
        array=~isnan(x_odd);
        edgeArray=diff([0; array; 0]);
        indices_odd=[find(edgeArray>0) find(edgeArray<0)-1];
        
        array=~isnan(x_even);
        edgeArray=diff([0;array;0]);
        indices_even=[find(edgeArray>0) find(edgeArray<0)-1];
        
        if ~isempty(indices_odd) && ~isempty(indices_even)
            for j=1:size(indices_odd,1)
                aa=indices_odd(j,1):indices_odd(j,2);
                for k=1:size(indices_even,1)
                    bb=indices_even(k,1):indices_even(k,2);
                    if ~((aa(end)<bb(1)) || (bb(end)<aa(1)))
                        XV=[x_odd(aa);flipud(x_even(bb));x_odd(aa(1))];
                        YV=[y_odd(aa); flipud(y_even(bb));y_odd(aa(1))];
                        XV_max=max(XV); XV_min=min(XV);
                        YV_max=max(YV); YV_min=min(YV);
                        [ptXtmp,ptYtmp]=meshgrid(XV_min:XV_max,YV_min:YV_max);
                        ptXtmp=ptXtmp(:);ptYtmp=ptYtmp(:);
                        if issorted(ptYtmp)
                            [in_on,~]=inpoly([ptXtmp ptYtmp],[XV YV]);
                        elseif issorted(ptXtmp)
                            [in_on,~]=inpoly([ptYtmp ptXtmp],[YV XV]);
                        else
                            [in_on,~]=inpoly([ptXtmp ptYtmp],[XV YV]);
                        end
                        pix_in_on=sub2ind([480,640],[ptYtmp(in_on==1);YV],[ptXtmp(in_on==1);XV]);
                        img_right(pix_in_on)=1;
                        break;
                    end
                end
            end
        end
    end
end