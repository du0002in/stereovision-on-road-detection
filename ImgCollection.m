terminate_all_threads;
folder_ind=1;
folder_exist=exist(['\Stereo\Image Collection(new Cam)' int2str(folder_ind)],'dir');
while folder_exist~=0
  folder_ind=folder_ind+1;
  folder_exist=exist(['\Stereo\Image Collection(new Cam)' int2str(folder_ind)],'dir');
end
mkdir(['\Stereo\Image Collection(new Cam)' int2str(folder_ind)]);
save_dir=['\Stereo\Image Collection(new Cam)' int2str(folder_ind)];

% % runningstatus(0);
% % closecommunication(0);
% % system('D:\Stereo\RealRun_Laptop_new_Cam\RealRun_Laptop_new_Cam_RS232\RealRun_Laptop_new_Cam_RS232\Debug\RealRun_Laptop_new_Cam_RS232.exe &');

[status, vid1, vid2]=openstereo;
if status==0
    disp ('Not able to initilize stereo cameras');
    return;
end
disp 'Stereo pair initiates sucessfully'
ii=0;
while stop_all_flag==2
    stop_all_flag=get_stop_all_flag;
    st=tic;
    ii=ii+1;
    trigger([vid1 vid2]);
    log_trigger_t(ii)=toc(tic);
    leftimg=getdata(vid1);
    rightimg=getdata(vid2);
    log_get_t(ii)=toc(tic);
    GA_inner_loop_logging('Ini');
    [act_v(ii) act_phai(ii)]=get_measurements;
    imwrite(leftimg,strcat(save_dir,'\left',int2str(ii),'.jpg'),'jpg');
    imwrite(rightimg,strcat(save_dir,'\right',int2str(ii),'.jpg'),'jpg');
    log_save_t(ii)=toc(st);
    pause(0.27);
    GA_inner_loop_logging('End');
    [loop_t, loop_phi, loop_v]=getloop_t_phi_v;
    logloop_data{ii}={loop_t,loop_phi,loop_v};
    ft=toc(st);
    act_t(ii)=ft;
end
terminate_all_threads;
stop([vid1 vid2]);
stoppreview([vid1 vid2]);
closepreview([vid1 vid2]);
A0_data_saving_date=datestr(now);
save([save_dir '\odemotry.mat']);






    
    
    
    
    
    
    