runningstatus(0);
closecommunication(0);
system('D:\Stereo\RealRun_Laptop_new_Cam\RealRun_Laptop_new_Cam_RS232\RealRun_Laptop_new_Cam_RS232\Debug\RealRun_Laptop_new_Cam_RS232.exe &');
disp 'Openning stereo pair'
[status, vid1, vid2]=openstereo;
if status==0
    disp ('Not able to initilize stereo cameras');
    return;
end
disp 'Stereo pair initiates sucessfully'

while (1)
    if get_runningstatus==0
    else
        status=run_RealRun_Laptop_new_Cam;
        if status~=1
            break;
        end
    end
    pause(0.05);
end

