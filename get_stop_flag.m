%flag=1 ==> keep on running
%flag=0 ==> stop the process
function flag=get_stop_flag
ele_no=0;
while(ele_no==0) %in case when the file is being written, it will read in nothing
    fid=fopen('D:\Stereo\GUI\runningflag.txt','rt');
    C=textscan(fid, '%d', 1);
    ele_no=numel(C{1});
    fclose(fid);
end
flag=C{1}(1);