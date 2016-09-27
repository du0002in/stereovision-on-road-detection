function [linelz, linerz, paraz, y_int]=P_Zebra_RANSAC_new(Datapointl, SubDatapointl, Datapointr, SubDatapointr, D_pre, H_pre, num_ite, focal, dist_left_right, cc_y)

[linelz, linerz, paraz]=P_Zebra_RANSAC_single_line(Datapointl, SubDatapointl, Datapointr, SubDatapointr, D_pre, H_pre, num_ite, focal, dist_left_right, cc_y);
if sum(linelz==1)<20 || sum(linerz==1)<20
    paraz=zeros(5,2);
    linelz=0; linerz=0;
    y_int=0;
    return;
end

cluster=C_DBScan1(zeros(1,sum(linelz==1)),Datapointl(1,linelz==1),15.0,2);
if max(cluster)==0
    paraz=zeros(5,2);
    linelz=0; linerz=0;
    y_int=0;
    return;
end
for i=1:max(cluster)
    log_no_in_cluster(i)=sum(cluster==i);
end
[~,i]=max(log_no_in_cluster);
indtmp=(cluster~=i);
indtmp2=find(linelz==1);
linelz(indtmp2(indtmp))=-1;

cluster=C_DBScan1(zeros(1,sum(linerz==1)),Datapointr(1,linerz==1),15.0,2);
if max(cluster)==0
    paraz=zeros(5,2);
    linelz=0; linerz=0;
    y_int=0;
    return;
end
for i=1:max(cluster)
    log_no_in_cluster(i)=sum(cluster==i);
end
[~,i]=max(log_no_in_cluster);
indtmp=(cluster~=i);
indtmp2=find(linerz==1);
linerz(indtmp2(indtmp))=-1;


if (sum(linelz==0)<10 || sum(linerz==0)<10) || sum(linelz==1)<20 || sum(linerz==1)<20
    paraz=zeros(5,2);
    linelz=0; linerz=0;
    y_int=0;
    return;
end

phi=atan(D_pre/focal);
theta_l=atan((-sin(phi)*paraz(2,1)+paraz(3,1)/2/focal*cos(phi))/(sin(phi)*sin(phi)-cos(phi)*cos(phi)));
theta_r=atan((-sin(phi)*paraz(4,1)+paraz(5,1)/2/focal*cos(phi))/(sin(phi)*sin(phi)-cos(phi)*cos(phi)));
max_y_l=max(Datapointl(1,linelz==1));
max_y_r=max(Datapointr(1,linerz==1));
max_y=min([max_y_l max_y_r]);
ind_l1=(Datapointl(1,:)>max_y_l) & (linelz==0);
ind_r1=(Datapointr(1,:)>max_y_r) & (linerz==0);

if (sum(ind_l1)<10 || sum(ind_r1)<10)
    paraz1=zeros(5,2); linelz1=0; linerz1=0; valid1=0;
else
    Dl1=Datapointl(:,ind_l1);
    Dr1=Datapointr(:,ind_r1);
    min_y_l_remain=min(Dl1(1,:));
    min_y_r_remain=min(Dr1(1,:));
    min_y_remain=max([min_y_l_remain min_y_r_remain]);
    delta=min_y_remain-max_y;
    sigma=max([delta/2 5]);
    miu=mean([min_y_remain max_y]);
    y_rnd=normrnd(miu,sigma,1,num_ite);
    x_rnd_l=-paraz(2,1)*y_rnd-paraz(3,1)/2;
    x_rnd_r=-paraz(4,1)*y_rnd-paraz(5,1)/2;
    Rnd_l=[y_rnd; x_rnd_l];
    Rnd_r=[y_rnd; x_rnd_r]; 
    [linelz1, linerz1, paraz1, break_loc1]=P_Zebra_RANSAC_2nd_line(Dl1,Dl1,Dr1,Dr1,D_pre, H_pre, num_ite, focal, dist_left_right, cc_y,Rnd_l, Rnd_r, theta_l, theta_r);
    if (sum(linelz1)>10 && sum(linerz1)>10)
        valid1=1;
    else valid1=0;
    end
end

min_y_l=min(Datapointl(1,linelz==1));
min_y_r=min(Datapointr(1,linerz==1));
min_y=max([min_y_l min_y_r]);
ind_l2=(Datapointl(1,:)<min_y_l) & (linelz==0);
ind_r2=(Datapointr(1,:)<min_y_r) & (linerz==0);
if (sum(ind_l2)<10 || sum(ind_r2)<10)
    paraz2=zeros(5,2); linelz2=0; linerz2=0; valid2=0;
else
    Dl2=Datapointl(:,ind_l2);
    Dr2=Datapointr(:,ind_r2);
    max_y_l_remain=max(Dl2(1,:));
    max_y_r_remain=max(Dr2(1,:));
    max_y_remain=min([max_y_l_remain max_y_r_remain]);
    delta=-max_y_remain+min_y;
    sigma=max([delta/2 5]);
    miu=mean([max_y_remain min_y]);
    y_rnd=normrnd(miu,sigma,1,num_ite);
    x_rnd_l=-paraz(2,1)*y_rnd-paraz(3,1)/2;
    x_rnd_r=-paraz(4,1)*y_rnd-paraz(5,1)/2;
    Rnd_l=[y_rnd; x_rnd_l];
    Rnd_r=[y_rnd; x_rnd_r];  
    [linelz2, linerz2, paraz2, break_loc2]=P_Zebra_RANSAC_2nd_line(Dl2,Dl2,Dr2,Dr2,D_pre, H_pre, num_ite, focal, dist_left_right, cc_y,Rnd_l, Rnd_r, theta_l, theta_r);
    if (sum(linelz2)>10 && sum(linerz2)>10)
        valid2=1;
    else valid2=0;
    end
end

if (valid1==0) && (valid2==0)
    paraz=zeros(5,2);
    linelz=0; linerz=0;
    y_int=0;
    return;
elseif ((valid1==1) && (valid2==0)) || ...
        ((valid1==1) && (valid2==1) && ((sum(linelz1)+sum(linerz1))>(sum(linelz2)+sum(linerz2))))
    paraz=[paraz paraz1];
    y_int=[break_loc1 break_loc1];
    indtmp= (linelz1==1);
    indtmp2=find(ind_l1==1);
    indtmp3=indtmp2(indtmp);
    linelz(indtmp3)=1;
    indtmp= (linerz1==1);
    indtmp2=find(ind_r1==1);
    indtmp3=indtmp2(indtmp);
    linerz(indtmp3)=1;
elseif ((valid1==0) && (valid2==1)) || ...
        ((valid1==1) && (valid2==1) && ((sum(linelz1)+sum(linerz1))<=(sum(linelz2)+sum(linerz2))))
    paraz=[paraz2 paraz];
    y_int=[break_loc2 break_loc2];
    indtmp= (linelz2==1);
    indtmp2=find(ind_l2==1);
    indtmp3=indtmp2(indtmp);
    linelz(indtmp3)=1;
    indtmp= (linerz2==1);
    indtmp2=find(ind_r2==1);
    indtmp3=indtmp2(indtmp);
    linerz(indtmp3)=1;
end
linelz(linelz==-1)=0;
linerz(linerz==-1)=0;
end