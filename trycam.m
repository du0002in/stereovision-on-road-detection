vid1=videoinput('winvideo', 1, 'RGB24_640x360');
    set(vid1,'FramesPerTrigger',1);
    set(vid1,'TriggerRepeat',Inf);
    triggerconfig(vid1,'manual');
    start(vid1);
preview(vid1);
for i=1:10
    trigger(vid1);
    leftimg=getdata(vid1);
    leftimg=double(max(leftimg,[],3));
    leftimg(leftimg>200)=255;
    leftimg(leftimg<=200)=0;
    imshow(leftimg,[]);
    pause(0.01);
end
stop(vid1);