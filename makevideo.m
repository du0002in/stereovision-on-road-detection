%making video
clear all;
f='D:\StereoLaptop_new\Image Collection(new Cam)91';
M = avifile([f,'\AfterRain','.avi'],'fps',6,'compression','None');%To build the video, with name and how many figures do u want in one second
% Figure here
log_ite_s=[2 461 1110 1490];
log_ite_e=[378 812 1444 1540];
%for Image Collection(new Cam)61 at dust
log_ite_s=[293];
log_ite_e=[1152];
%for Image Collection(new Cam)90 at night
log_ite_s=[216];
log_ite_e=[942];
%for Image Collection(new Cam)91 after rain
log_ite_s=[282];
log_ite_e=[993];

for outloop=1:1
ite_s=log_ite_s(outloop);
ite_e=log_ite_e(outloop);
for i=ite_s:ite_e
    img=strcat(f,'\3DFresult_letter',int2str(i),'.jpg');
    if exist(img)
    imshow(img);
    frame = getframe(gcf);
    M = addframe(M,frame);
    end
end
end
M=close(M);