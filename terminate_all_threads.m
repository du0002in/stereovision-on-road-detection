i=1;
try
    GA_inner_loop_logging('End');
    X=[int2str(i), '. GA_inner_loop_logging terminated'];
    disp(X);
catch
    X=[int2str(i), '. GA_inner_loop_logging not initilized'];
    disp(X);
end

i=i+1;
try
    get_status('End');
    X=[int2str(i), '. get_status terminated'];
    disp(X);
catch
    X=[int2str(i), '. get_status not initilized'];
    disp(X);
end

i=i+1;
try
    RealRunPXI_inner_loop_logging_w_initial('End',[0,0,0]);
    X=[int2str(i), '. RealRunPXI_inner_loop_logging_w_initial terminated'];
    disp(X);
catch
    X=[int2str(i), '. RealRunPXI_inner_loop_logging_w_initial not initilized'];
    disp(X);
end

i=i+1;
try
    C_MPCsendtarget('End',zeros(2,15),0.1);
    X=[int2str(i), '. C_MPCsendtarget terminated'];
    disp(X);
catch
    X=[int2str(i), '. C_MPCsendtarget not initilized'];
    disp(X);
end

i=i+1;
try
    GA_inner_loop_logging_w_initial('End',[0,0,0]);
    X=[int2str(i), '. GA_inner_loop_logging_w_initial terminated'];
    disp(X);
catch
    X=[int2str(i), '. GA_inner_loop_logging_w_initial not initilized'];
    disp(X);
end

i=i+1;
try
    GA_VehicleSimulatorMPC('End',0, 0, 0, 0, 0, 0, 0, [0 0 0]);
    X=[int2str(i), '. GA_VehicleSimulatorMPC terminated'];
    disp(X);
catch
    X=[int2str(i), '. GA_VehicleSimulatorMPC not initilized'];
    disp(X);
end