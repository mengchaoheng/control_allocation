clear all;
close all;
clc;
clear all;
close all;
clc;
addpath(genpath(pwd));
% you can run on terminal 
% ulog2csv log_8_2021-5-20-11-52-08.ulg 
% to get csv files
% =====================1==========================
% Install pyulog using pip first.https://github.com/PX4/pyulog.
% in MacOS, it maybe have been installed by the px4-dev
% =====================2==========================
% Make sure it has installed ulog2csv correctly (check the output of which ulog2csv in Linux/MacOS or where ulog2csv in Windows).
% =====================3==========================
% Change the following line in ulogviewver.m:
% command = ['!/usr/local/bin/ulog2csv ' ulgFileName '.ulg'];
% to 
% command = ['!your ulog2csv path' ulgFileName '.ulg'];
% and 
% ulgFileName = '00_41_22'; 
% to 
% ulgFileName = 'your log name'; 

% ----fig size, you have to change it for your fig

d2r=pi/180;
r2d=180/pi;
%% new method: defualt pid for endurance load wp , remove I for wind and max vel
% endurance:(17_14_37 1127s) 12_07_46 1137s   11_10_22 1140s  11_37_27 1120s  
% load:06_28_53 06_45_38  06_54_35
% wp:17_09_30 
% wind: 11_39_14  15_57_58
 % mav vel: 09_19_16 10_27_37

%% 2023 flight data
% Users/mch/Documents/FlightLog/log/2023-02-24/06_12_56.ulg
% 06_22_46
% 06_44_14
% 07_04_47
% 07_26_57

% 03_43_50
% /Users/mch/Documents/FlightLog/log/2023-02-23/07_44_32.ulg
% 07_45_53
% 07_48_07
% 07_50_17
% 08_06_35
% 08_29_05
% 08_39_27

% 20240614 08_15_39 08_39_28 09_05_41 09_15_41 
% 20240616  05_31_54  14_03_42 06_20_57 15_12_39
 %% pca test:16_37_44 12_15_41 15_03_31  df-1.12.3 (15_10_14)

 % real flight 07_18_37  12_03_09   (important high rate:08_43_49  08_44_10 all:08_21_58 ekf:08_49_36)
ulgFileName = '11_22_05'; % the ulog file name. load 06_28_53 06_45_38  06_54_35 wind 11_39_14  15_57_58 wp:17_09_30, endurance:17_14_37 (18.9) 06_26_32 (17.8) 06_50_20 (19.6) (19.5) 07_10_30 (18.6)
tmp=[ ulgFileName '.mat'];
% exist tmp var05_31_54.ulg
if exist(tmp,"file")
    load(ulgFileName,'log');
else
    if ismac
        % on macOS, run " which ulog2csv " on terminal to get it.
        command = ['!/Users/mch/opt/anaconda3/bin/ulog2csv ' ulgFileName '.ulg']; % /usr/local/bin/ is the path of ulog2csv, 
    else
        % on windows and linux just make sure you have installed pyulog
        command = ['!ulog2csv ' ulgFileName '.ulg']; % have installed ulog2csv,
    end

	eval(command);
    log.data = csv_topics_to_d(ulgFileName);
    log.FileName = ulgFileName;
    log.version = 1.0;
    log.params = '';
    log.messages = '';
    log.info = '';
    %run add_fields_in_preprocessing.m
    save(ulgFileName,'log')
    delete(['*' ulgFileName '*.csv'])
end
%%
if(isfield(log.data, 'vehicle_angular_velocity_0'))
    vehicle_angular_velocity=log.data.vehicle_angular_velocity_0{:,:};
    [rate_N,~]=size(vehicle_angular_velocity(:,1));
    rate_delta_t=zeros(rate_N-1,1);
    for i=1:rate_N-1
        rate_delta_t(i)=(vehicle_angular_velocity(i+1,1))*1e-6-(vehicle_angular_velocity(i,1))*1e-6;
    end
end 
if(isfield(log.data, 'vehicle_angular_acceleration_0'))
    vehicle_angular_acceleration=log.data.vehicle_angular_acceleration_0{:,:};
    [rate_acc_N,~]=size(vehicle_angular_acceleration(:,1));
    rate_acc_delta_t=zeros(rate_acc_N-1,1);
    for i=1:rate_acc_N-1
        rate_acc_delta_t(i)=(vehicle_angular_acceleration(i+1,1))*1e-6-(vehicle_angular_acceleration(i,1))*1e-6;
    end
end 
if(isfield(log.data, 'vehicle_rates_setpoint_0'))
    vehicle_rates_setpoint=log.data.vehicle_rates_setpoint_0{:,:};
    % [rate_N,~]=size(vehicle_rates_setpoint(:,1));
    % rate_delta_t=zeros(rate_N-1,1);
    % for i=1:rate_N-1
    %     rate_delta_t(i)=(vehicle_rates_setpoint(i+1,1))*1e-6-(vehicle_rates_setpoint(i,1))*1e-6;
    % end
end 
if(isfield(log.data, 'vehicle_attitude_0'))
    vehicle_attitude=log.data.vehicle_attitude_0{:,:};
    q_0=vehicle_attitude(:,3);
    q_1=vehicle_attitude(:,4);
    q_2=vehicle_attitude(:,5);
    q_3=vehicle_attitude(:,6);
    Roll=quat_to_roll([q_0 q_1 q_2 q_3]);
    Pitch=quat_to_pitch([q_0 q_1 q_2 q_3]);
    Yaw=quat_to_yaw([q_0 q_1 q_2 q_3]);
end 
if(isfield(log.data, 'vehicle_attitude_setpoint_0'))
    vehicle_attitude_setpoint=log.data.vehicle_attitude_setpoint_0{:,:};
    q_0_setpoint=vehicle_attitude_setpoint(:,6);
    q_1_setpoint=vehicle_attitude_setpoint(:,7);
    q_2_setpoint=vehicle_attitude_setpoint(:,8);
    q_3_setpoint=vehicle_attitude_setpoint(:,9);
    Roll_setpoint=quat_to_roll([q_0_setpoint q_1_setpoint q_2_setpoint q_3_setpoint]);
    Pitch_setpoint=quat_to_pitch([q_0_setpoint q_1_setpoint q_2_setpoint q_3_setpoint]);
    Yaw_setpoint=quat_to_yaw([q_0_setpoint q_1_setpoint q_2_setpoint q_3_setpoint]);
    [attitude_N,~]=size(vehicle_attitude_setpoint(:,1));
    attitude_delta_t=zeros(attitude_N-1,1);
    for i=1:attitude_N-1
        attitude_delta_t(i)=(vehicle_attitude_setpoint(i+1,1))*1e-6-(vehicle_attitude_setpoint(i,1))*1e-6;
    end
end 

if(isfield(log.data, 'vehicle_local_position_0'))
    vehicle_local_position=log.data.vehicle_local_position_0{:,:};
    XYZ=vehicle_local_position(:,6:8);
    V_XYZ=vehicle_local_position(:,12:14);
end 

if(isfield(log.data, 'vehicle_local_position_setpoint_0'))
    vehicle_local_position_setpoint=log.data.vehicle_local_position_setpoint_0{:,:};
    XYZ_setpoint=vehicle_local_position_setpoint(:,2:4);
    V_XYZ_setpoint=vehicle_local_position_setpoint(:,7:9);
    [pose_N,~]=size(vehicle_local_position_setpoint(:,1));
    pose_delta_t=zeros(pose_N-1,1);
    for i=1:pose_N-1
        pose_delta_t(i)=(vehicle_local_position_setpoint(i+1,1))*1e-6-(vehicle_local_position_setpoint(i,1))*1e-6;
    end
end 




if(isfield(log.data, 'actuator_controls_0_0'))
    actuator_controls=log.data.actuator_controls_0_0{:,:};   
    Roll_control=actuator_controls(:,3);
    Pitch_control=actuator_controls(:,4);
    Yaw_control=actuator_controls(:,5);
    if(ismember('indi_fb_0_', log.data.actuator_controls_0_0.Properties.VariableNames))
        indi_feedback=actuator_controls(:,11:13);
        error_feedback=actuator_controls(:,14:16);
    end
    [actuator_N,~]=size(actuator_controls(:,1));
    actuator_delta_t=zeros(actuator_N-1,1);
    for i=1:actuator_N-1
        actuator_delta_t(i)=(actuator_controls(i+1,1))*1e-6-(actuator_controls(i,1))*1e-6;
    end
end
if(isfield(log.data, 'actuator_outputs_0'))
    actuator_outputs=log.data.actuator_outputs_0{:,:};   
    cs1=actuator_outputs(:,7);
    cs2=actuator_outputs(:,8);
    cs3=actuator_outputs(:,9);
    cs4=actuator_outputs(:,10);
    [cs_N,~]=size(actuator_outputs(:,1));
    cs_delta_t=zeros(cs_N-1,1);
    for i=1:cs_N-1
        cs_delta_t(i)=(actuator_outputs(i+1,1))*1e-6-(actuator_outputs(i,1))*1e-6;
    end
end


if(isfield(log.data, 'allocation_value_0'))
    allocation_value=log.data.allocation_value_0{:,:};

end 

if(isfield(log.data, 'vehicle_visual_odometry_0'))
    vehicle_visual_odometry=log.data.vehicle_visual_odometry_0{:,:};
    visual_odometry_X=vehicle_visual_odometry(:,3);
    visual_odometry_Y=vehicle_visual_odometry(:,4);
    visual_odometry_Z=vehicle_visual_odometry(:,5);
    visual_odometry_q0=vehicle_visual_odometry(:,6);
    visual_odometry_q1=vehicle_visual_odometry(:,7);
    visual_odometry_q2=vehicle_visual_odometry(:,8);
    visual_odometry_q3=vehicle_visual_odometry(:,9);
    
end


%% 
if(isfield(log.data, 'vehicle_angular_velocity_0') && isfield(log.data, 'vehicle_rates_setpoint_0'))
    fig1=figure(1);
    subplot(311)
    plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,2)*r2d,'k-','LineWidth',1);hold on;
    plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,3)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -100 100]);
    xlabel({'时间(s)'});
    ylabel('滚转(deg/s)')
    title('角速度');
    legend('给定','响应');
    %% 
    subplot(312)
    plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,3)*r2d,'k-','LineWidth',1);hold on;
    plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,4)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -100 100]);
    xlabel({'时间(s)'});
    ylabel('俯仰(deg/s)')
    % title('Pitch angular velocity');
    legend('给定','响应');
    %% 
    subplot(313)
    plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,4)*r2d,'k-','LineWidth',1);hold on;
    plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,5)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -100 100]);
    xlabel({'时间(s)'});
    ylabel('偏航(deg/s)')
    % title('Yaw angular velocity');
    legend('给定','响应');
    %% 
    PlotToFileColorPDF(fig1,'results/pqr.png',15,20); % or 'pqr.png'
end


%% and maybe more figure, all in the variable "log.data"
if(isfield(log.data, 'vehicle_attitude_setpoint_0') && isfield(log.data, 'vehicle_attitude_0'))
    fig2=figure(2);
    subplot(311)
    plot((vehicle_attitude_setpoint(:,1))*1e-6, Roll_setpoint*r2d,'k-','LineWidth',1);hold on;
    plot((vehicle_attitude(:,1))*1e-6, Roll*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -100 100]);
    xlabel({'时间(s)'});
    ylabel('滚转角(deg)')
    title('姿态角');
    legend('给定','响应');
    %% and maybe more figure, all in the variable "log.data"
    subplot(312)
    plot((vehicle_attitude_setpoint(:,1))*1e-6, Pitch_setpoint*r2d,'k-','LineWidth',1);hold on;
    plot((vehicle_attitude(:,1))*1e-6, Pitch*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -100 100]);
    xlabel({'时间(s)'});
    ylabel('俯仰角(deg)')
    legend('给定','响应');
    %% and maybe more figure, all in the variable "log.data"
    subplot(313)
    plot((vehicle_attitude_setpoint(:,1))*1e-6, Yaw_setpoint*r2d,'k-','LineWidth',1);hold on;
    plot((vehicle_attitude(:,1))*1e-6, Yaw*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -100 100]);
    xlabel({'时间(s)'});
    ylabel('偏航角(deg)')
    legend('给定','响应');
    %% 
    PlotToFileColorPDF(fig2,'results/RPY.png',15,20); % or 'RPY.png'
end





%% and maybe more figure, all in the variable "log.data"
if(isfield(log.data, 'vehicle_local_position_0') && isfield(log.data, 'vehicle_local_position_setpoint_0'))

    fig3=figure(3);
    subplot(311)
    plot((vehicle_local_position_setpoint(:,1))*1e-6, V_XYZ_setpoint(:,1),'k-','LineWidth',1);hold on;
    plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,1),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    title('速度');
    xlabel({'时间(s)'});
    ylabel('V_X(m/s)')
    legend('给定','响应');
    %% and maybe more figure, all in the variable "log.data"
    subplot(312)
    plot((vehicle_local_position_setpoint(:,1))*1e-6, V_XYZ_setpoint(:,2),'k-','LineWidth',1);hold on;
    plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,2),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('V_Y(m/s)')
    legend('给定','响应');
    %% and maybe more figure, all in the variable "log.data"
    subplot(313)
    plot((vehicle_local_position_setpoint(:,1))*1e-6, V_XYZ_setpoint(:,3),'k-','LineWidth',1);hold on;
    plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,3),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('V_Z(m/s)')
    legend('给定','响应');
    PlotToFileColorPDF(fig3,'results/V_XYZ.png',15,20);
elseif(isfield(log.data, 'vehicle_local_position_0'))
    fig3=figure(3);
    subplot(311)
    plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,1),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    title('速度');
    xlabel({'时间(s)'});
    ylabel('V_X(m/s)')
    %% and maybe more figure, all in the variable "log.data"
    subplot(312)
    plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,2),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('V_Y(m/s)')
    %% and maybe more figure, all in the variable "log.data"
    subplot(313)
    plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,3),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('V_Z(m/s)')
    PlotToFileColorPDF(fig3,'results/V_XYZ.png',15,20);
end
%% 


if(isfield(log.data, 'vehicle_local_position_0') && isfield(log.data, 'vehicle_local_position_setpoint_0'))

%% and maybe more figure, all in the variable "log.data"
    fig4=figure(4);
    subplot(311)
    plot((vehicle_local_position_setpoint(:,1))*1e-6, XYZ_setpoint(:,1),'k-','LineWidth',1);hold on;
    plot((vehicle_local_position(:,1))*1e-6, XYZ(:,1),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    title('位置');
    xlabel({'时间(s)'});
    ylabel('X(m)')
    legend('给定','响应');
    %% and maybe more figure, all in the variable "log.data"
    subplot(312)
    plot((vehicle_local_position_setpoint(:,1))*1e-6, XYZ_setpoint(:,2),'k-','LineWidth',1);hold on;
    plot((vehicle_local_position(:,1))*1e-6, XYZ(:,2),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('Y(m)')
    legend('给定','响应');
    %% and maybe more figure, all in the variable "log.data"
    subplot(313)
    plot((vehicle_local_position_setpoint(:,1))*1e-6, XYZ_setpoint(:,3),'k-','LineWidth',1);hold on;
    plot((vehicle_local_position(:,1))*1e-6, XYZ(:,3),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('Z(m)')
    legend('给定','响应');
    PlotToFileColorPDF(fig4,'results/P_XYZ.png',15,20);
elseif(isfield(log.data, 'vehicle_local_position_0'))
    fig4=figure(4);
    subplot(311)
    plot((vehicle_local_position(:,1))*1e-6, XYZ(:,1),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    title('位置');
    xlabel({'时间(s)'});
    ylabel('X(m)')
    %% and maybe more figure, all in the variable "log.data"
    subplot(312)
    plot((vehicle_local_position(:,1))*1e-6, XYZ(:,2),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('Y(m)')
    %% and maybe more figure, all in the variable "log.data"
    subplot(313)
    plot((vehicle_local_position(:,1))*1e-6, XYZ(:,3),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('Z(m)')
    PlotToFileColorPDF(fig4,'results/P_XYZ.png',15,20);
%%
end
%% 
% 



%% 
if(isfield(log.data, 'actuator_controls_0_0'))

    fig5=figure(5);
    subplot(411)
    plot((actuator_controls(:,1))*1e-6, Roll_control(:,1),'r-','LineWidth',1);hold on;
    plot((actuator_controls(:,1))*1e-6, Pitch_control(:,1),'k--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    plot((actuator_controls(:,1))*1e-6, Yaw_control(:,1),'b-.','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    title('控制器输出');
    xlabel({'时间(s)'});
    ylabel('控制量')
    legend('滚转','俯仰',  '偏航');
    subplot(412)
    plot((actuator_controls(:,1))*1e-6, Roll_control(:,1),'r-','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('滚转')
    subplot(413)
    plot((actuator_controls(:,1))*1e-6, Pitch_control(:,1),'k--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('俯仰')
    subplot(414)
    plot((actuator_controls(:,1))*1e-6, Yaw_control(:,1),'b-.','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('偏航')
    %% 
    PlotToFileColorPDF(fig5,'results/actuator_controls.png',15,20);


%% 
    if(ismember('indi_fb_0_', log.data.actuator_controls_0_0.Properties.VariableNames))
    fig6=figure(6);
    subplot(411)
    plot((actuator_controls(:,1))*1e-6, indi_feedback(:,1),'r-','LineWidth',1);hold on;
    plot((actuator_controls(:,1))*1e-6, indi_feedback(:,2),'k--','LineWidth',1);hold on;
    plot((actuator_controls(:,1))*1e-6, indi_feedback(:,3),'b-.','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    title('indi控制量');
    xlabel({'时间(s)'});
    ylabel('控制量')
    legend('滚转','俯仰',  '偏航');
    subplot(412)
    plot((actuator_controls(:,1))*1e-6, indi_feedback(:,1),'r-','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('滚转')
    subplot(413)
    plot((actuator_controls(:,1))*1e-6, indi_feedback(:,2),'k--','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('俯仰')
    subplot(414)
    plot((actuator_controls(:,1))*1e-6, indi_feedback(:,3),'b-.','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('偏航')
    % %% 
    PlotToFileColorPDF(fig6,'results/indi_feedback.png',15,20);


    fig7=figure(7);
    subplot(411)
    plot((actuator_controls(:,1))*1e-6, error_feedback(:,1),'r-','LineWidth',1);hold on;
    plot((actuator_controls(:,1))*1e-6, error_feedback(:,2),'k--','LineWidth',1);hold on;
    plot((actuator_controls(:,1))*1e-6, error_feedback(:,3),'b-.','LineWidth',1);hold on;

    grid on;
    % axis([-inf inf -0.5 0.5]);
    title('误差反馈控制量');
    xlabel({'时间(s)'});
    ylabel('控制量')
    legend('滚转','俯仰',  '偏航');
    subplot(412)
    plot((actuator_controls(:,1))*1e-6, error_feedback(:,1),'r-','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('滚转')
    subplot(413)
    plot((actuator_controls(:,1))*1e-6, error_feedback(:,2),'k--','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('俯仰')
    subplot(414)
    plot((actuator_controls(:,1))*1e-6, error_feedback(:,3),'b-.','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf -0.5 0.5]);
    xlabel({'时间(s)'});
    ylabel('偏航')
    % %% 
    PlotToFileColorPDF(fig7,'results/error_control_input.png',15,20);
    end

end




%% 
if(isfield(log.data, 'actuator_outputs_0'))

    fig8=figure(8);
    subplot(511)
    plot((actuator_outputs(:,1))*1e-6, cs1(:,1),'r-','LineWidth',1);hold on;
    plot((actuator_outputs(:,1))*1e-6, cs2(:,1),'k--','LineWidth',1);hold on;
    plot((actuator_outputs(:,1))*1e-6, cs3(:,1),'b-.','LineWidth',1);hold on;
    plot((actuator_outputs(:,1))*1e-6, cs4(:,1),'g-','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf 1000 2000]);
    title('操作面指令');
    xlabel({'时间(s)'});
    ylabel('偏转指令(pwm)')
    legend('操作面1','操作面2','操作面3','操作面4');
    subplot(512)
    plot((actuator_outputs(:,1))*1e-6, cs1(:,1),'r-','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf 1000 2000]);
    xlabel({'时间(s)'});
    ylabel('操作面1')
    subplot(513)
    plot((actuator_outputs(:,1))*1e-6, cs2(:,1),'k--','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf 1000 2000]);
    xlabel({'时间(s)'});
    ylabel('操作面2')
    subplot(514)
    plot((actuator_outputs(:,1))*1e-6, cs3(:,1),'b-.','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf 1000 2000]);
    xlabel({'时间(s)'});
    ylabel('操作面3')
    subplot(515)
    plot((actuator_outputs(:,1))*1e-6, cs4(:,1),'g-','LineWidth',1);hold on;
    grid on;
    % axis([-inf inf 1000 2000]);
    xlabel({'时间(s)'});
    ylabel('操作面4')
    %% 
    PlotToFileColorPDF(fig8,'results/cs.png',15,20);
end




%%
if(isfield(log.data, 'vehicle_local_position_0') && isfield(log.data, 'vehicle_local_position_setpoint_0'))

    % 创建图形窗口
    fig9=figure(9);

    % 绘制第一条轨迹
    plot3(XYZ_setpoint(:,1), XYZ_setpoint(:,2), -XYZ_setpoint(:,3), 'LineStyle', '-', 'LineWidth', 1);
    hold on;
    % 绘制第二条轨迹
    plot3(XYZ(:,1), XYZ(:,2), -XYZ(:,3), 'LineStyle', ':', 'LineWidth', 1);

    % 设置标题和轴标签
    title('轨迹跟踪');
    xlabel('X(m)');
    ylabel('Y(m)');
    zlabel('Z(m)');
    axis([-7 7 -7 7 -0.5 5]);% indoor
    % 设置图例
    legend1=legend('给定轨迹', '位置响应');
    set(legend1,...
        'Position',[0.154729910714286 0.304791666666668 0.146450892857143 0.0916071428571428]);
    % 设置网格
    grid on;

    % 设置视角
    view(45, 30);
    hold off;
    PlotToFileColorPDF(fig9,'results/trj.png',20,20);
end



%% 
if(isfield(log.data, 'vehicle_angular_velocity_0'))
    figure(10);
    plot(1:rate_N-1, rate_delta_t,'k-','LineWidth',1);hold on;
    ylabel('time (s)');
    title('vehicle angular velocity');
    disp('mean(rate_delta_t)');
    mean(rate_delta_t)
end
if(isfield(log.data, 'vehicle_attitude_setpoint_0'))
    figure(11);
    plot(1:attitude_N-1, attitude_delta_t,'k-','LineWidth',1);hold on;
    ylabel('time (s)');
    title('vehicle attitude setpoint');
    disp('mean(attitude_delta_t)');
    mean(attitude_delta_t)

end
if(isfield(log.data, 'vehicle_local_position_setpoint_0'))
    figure(12);
    plot(1:pose_N-1, pose_delta_t,'k-','LineWidth',1);hold on;
    ylabel('time (s)');
    title('vehicle local position setpoint');
    disp('mean(pose_delta_t)');
    mean(pose_delta_t)
end
if(isfield(log.data, 'actuator_controls_0_0'))
    figure(13);
    plot(1:actuator_N-1, actuator_delta_t,'k-','LineWidth',1);hold on;
    ylabel('time (s)');
    title('actuator controls');
    disp('mean(actuator_delta_t)');
    mean(actuator_delta_t)
end 
if(isfield(log.data, 'actuator_outputs_0'))
    figure(14);
    plot(1:cs_N-1, cs_delta_t,'k-','LineWidth',1);hold on;
    ylabel('time (s)');
    title('actuator outputs');
    disp('mean(cs_delta_t)');
    mean(cs_delta_t)
end

if(isfield(log.data, 'vehicle_angular_acceleration_0'))
    fig15=figure(15);
    subplot(311)
    plot((vehicle_angular_acceleration(:,1))*1e-6, vehicle_angular_acceleration(:,3),'k-','LineWidth',1);hold on;
    plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,3),'--','LineWidth',1,'color',[0.6,0.2,0,0.5]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    title('角加速度');
    xlabel({'时间(s)'});
    ylabel('p(rad/s^2)')
    legend('angular acc','gyro');
    %% and maybe more figure, all in the variable "log.data"
    subplot(312)
    plot((vehicle_angular_acceleration(:,1))*1e-6, vehicle_angular_acceleration(:,4),'k-','LineWidth',1);hold on;
    plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,4),'--','LineWidth',1,'color',[0.6,0.2,0,0.5]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('q(rad/s^2)')
    legend('angular acc','gyro');
    %% and maybe more figure, all in the variable "log.data"
    subplot(313)
    plot((vehicle_angular_acceleration(:,1))*1e-6, vehicle_angular_acceleration(:,5),'k-','LineWidth',1);hold on;
    plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,5),'--','LineWidth',1,'color',[0.6,0.2,0,0.5]);hold on;
    grid on;
    % axis([0 1200 -inf inf]);
    xlabel({'时间(s)'});
    ylabel('r(rad/s^2)')
    legend('angular acc','gyro');
    PlotToFileColorPDF(fig15,'results/vehicle_angular_acceleration.png',15,20);

end

if(isfield(log.data, 'vehicle_visual_odometry_0') && isfield(log.data, 'vehicle_attitude_0')) 
    fig16=figure(16);
    plot((vehicle_visual_odometry(:,1))*1e-6, visual_odometry_q0,'k-','LineWidth',1);hold on;
    plot((vehicle_attitude(:,1))*1e-6, q_0,'r-','LineWidth',1);hold on;grid on;
    % 采样频率不一样，只能粗略估计
    %计算时延方法1
    % % 使用 finddelay 函数
    % timeDelay = finddelay(q_0, visual_odometry_q0);
    % % 显示结果
    % disp(['The time delay between the signals is ', num2str(timeDelay), ' samples.']);
    %计算时延方法2
    % % 计算互相关
    % [c, lags] = xcorr(visual_odometry_q0, q_0);
    % 
    % % 找到最大相关性的位置
    % [~, I] = max(c);
    % timeDelay = lags(I);
    % 
    % % 显示结果
    % disp(['The time delay between the signals is ', num2str(timeDelay), ' samples.']);
    %计算时延方法3
    % X_q_0 = fft(q_0);
    % Y_visual_odometry_q0 = fft(visual_odometry_q0);
    % % 计算相位差
    % dPhi = angle(Y_visual_odometry_q0 ./ X_q_0);
    % % 计算时延
    % frequencies = (0:length(t)-1) * (fs/length(t));
    % timeDelay = mean(dPhi ./ (2 * pi * frequencies));
    % 
    % % 显示结果
    % disp(['The estimated time delay between the signals is ', num2str(timeDelay), ' seconds.']);

end




% [len1,~]=size(sensor_combined);
% [len2,~]=size(vehicle_visual_odometry);
% % [len2,~]=size(vehicle_gps_position);
% [len3,~]=size(vehicle_magnetometer);
% [len4,~]=size(vehicle_air_data);
% t1=1:1:len1;
% t2=1:1:len2;
% t3=1:1:len3;
% t4=1:1:len4;
% figure,
% plot(t1,sensor_combined(:,1),'k-');hold on;
% plot(t2,vehicle_visual_odometry(:,1),'r-');hold on;
% plot(t3,vehicle_magnetometer(:,1),'b');hold on;
% plot(t4,vehicle_air_data(:,1),'k--');hold on;
% legend('sensor combined','vehicle visual odometry','vehicle magnetometer','vehicle air data');
% figure,
% plot((vehicle_attitude(start:end,1))*1e-6,vehicle_attitude(start:end,4),'k:.');hold on;
% plot((vehicle_visual_odometry(start:end,1))*1e-6,vehicle_visual_odometry(start:end,7),'r:.');hold on;
% legend('vehicle attitude','vehicle visual odometry');



command_px4=log.data.actuator_controls_0_0{:,3:5}; % controller output
u_px4=log.data.allocation_value_0{:,12:15}; 

u_time_us=log.data.allocation_value_0{:,1}; 
u_delta_t=u_time_us(2:end)-u_time_us(1:end-1);
u_delta_t_s=u_delta_t*1e-6;


controls_time_us=log.data.actuator_controls_0_0{:,1}; 
controls_delta_t=controls_time_us(2:end)-controls_time_us(1:end-1);
controls_delta_t_s=controls_delta_t*1e-6;
% 指定保存路径及文件名
filename = 'flight.mat';

% 调用 save 函数进行保存
save(filename, 'command_px4','u_px4','controls_delta_t_s','u_delta_t_s');

del=log.data.allocation_value_0{:,1};
Delta_time=del*1e-6;mean(Delta_time)