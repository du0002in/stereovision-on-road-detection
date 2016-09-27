clear all;
close all;
runningstatus(1);
[status, vid1, vid2]=openstereo;
if status==0
    return;
end
fid=fopen('D:\Stereo\imglog.txt','w');
fclose(fid);
i=1;j=1;
while (get_runningstatus==1)
    t_start=tic;
    [v, phai]=get_measurements;
    
    trigger([vid1 vid2]);
    Img1=getdata(vid1);
    Img2=getdata(vid2);
%     logtmpg(i)=toc(t_start);
    imwrite(Img1,strcat('D:\Stereo\Image Collection(new Cam)12\left',int2str(i),'.jpg'),'jpg');
    imwrite(Img2,strcat('D:\Stereo\Image Collection(new Cam)12\right',int2str(i),'.jpg'),'jpg');
    t_end=toc(t_start);
    if t_end<0.36
        pause(0.36-t_end);
        t_end=toc(t_start);
    end
    fid=fopen('D:\Stereo\imglog.txt','a');
    formatSpec='%.3f   %4.3f   %4.3f\r\n';
    fprintf(fid,formatSpec,t_end, v, phai);
    fclose(fid);
    i=i+1;
end
stop([vid1 vid2]);