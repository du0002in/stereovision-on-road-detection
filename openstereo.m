function [status, varargout]=openstereo
imaqreset;
%left cam vid1
try
%     vid1=videoinput('winvideo', 3, 'BY8 _640x480','Tag','left_cam');
    vid1 = videoinput('tisimaq_r2013', 2, 'BY8 (640x480)', 'Tag', 'left_cam');
    set(vid1,'FramesPerTrigger',1);
    set(vid1,'TriggerRepeat',Inf);
    src1 = getselectedsource(vid1);
    src1.FrameRate='60.00';
    src1.ExposureAuto = 'Off';
    src1.ExposureAutoMaxValue = 1/30.0;
    triggerconfig(vid1,'manual');
catch err
    status=0;
    error(err);
    return;
end

%right cam vid2
try
%     vid2=videoinput('winvideo', 2, 'BY8 _640x480','Tag','right_cam');
    vid2 = videoinput('tisimaq_r2013', 1, 'BY8 (640x480)', 'Tag', 'right_cam');
    set(vid2,'FramesPerTrigger',1);
    set(vid2,'TriggerRepeat',Inf);
    src2 = getselectedsource(vid2);
    src2.FrameRate='60.00';
    src2.ExposureAuto = 'Off';
    src2.ExposureAutoMaxValue = 1/30.0;
    triggerconfig(vid2,'manual');
catch err
    stop(vid1);
    status=0;
    error(err);
    return;
end
status=1;
varargout{1}=vid1;
varargout{2}=vid2;
varargout{3}=src1;
varargout{4}=src2;
start([vid1 vid2]);
preview([vid1 vid2]);