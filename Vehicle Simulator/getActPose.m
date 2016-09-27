function [x, z, theta, phi, xc, d_theta, v]=getActPose
ele_no=0;
while(ele_no==0) %in case when the file is being written, it will read in nothing
    fid=fopen('D:\Stereo\Vehicle Simulator\act_x.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
x=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\Vehicle Simulator\act_z.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
z=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\Vehicle Simulator\act_theta.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
theta=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\Vehicle Simulator\act_phi.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
phi=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\Vehicle Simulator\act_xc.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
xc=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\Vehicle Simulator\act_d_theta.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
d_theta=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\Vehicle Simulator\act_v.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
v=C{1}(1);

end