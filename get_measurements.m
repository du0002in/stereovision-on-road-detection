function [v, phai]=get_measurements
ele_no=0;
while(ele_no==0) %in case when the file is being written, it will read in nothing
    fid=fopen('D:\Stereo\v_m.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
v=C{1}(1);

ele_no=0;
while(ele_no==0)
    fid=fopen('D:\Stereo\phai_m.txt','rt');
    C=textscan(fid, '%f', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
phai=C{1}(1);
end