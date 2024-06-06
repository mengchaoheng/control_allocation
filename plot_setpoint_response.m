% clear all;
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
%%
ulgFileName = '04_09_33'; % the ulog file name 
tmp=[ ulgFileName '.mat'];
% exist tmp var
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

vehicle_angular_velocity=log.data.vehicle_angular_velocity_0{:,:};
vehicle_rates_setpoint=log.data.vehicle_rates_setpoint_0{:,:};
vehicle_attitude=log.data.vehicle_attitude_0{:,:};
vehicle_attitude_setpoint=log.data.vehicle_attitude_setpoint_0{:,:};
q_0=vehicle_attitude(:,3);
q_1=vehicle_attitude(:,4);
q_2=vehicle_attitude(:,5);
q_3=vehicle_attitude(:,6);
Roll=quat_to_roll([q_0 q_1 q_2 q_3]);
Pitch=quat_to_pitch([q_0 q_1 q_2 q_3]);
Yaw=quat_to_yaw([q_0 q_1 q_2 q_3]);
q_0_setpoint=vehicle_attitude_setpoint(:,6);
q_1_setpoint=vehicle_attitude_setpoint(:,7);
q_2_setpoint=vehicle_attitude_setpoint(:,8);
q_3_setpoint=vehicle_attitude_setpoint(:,9);
Roll_setpoint=quat_to_roll([q_0_setpoint q_1_setpoint q_2_setpoint q_3_setpoint]);
Pitch_setpoint=quat_to_pitch([q_0_setpoint q_1_setpoint q_2_setpoint q_3_setpoint]);
Yaw_setpoint=quat_to_yaw([q_0_setpoint q_1_setpoint q_2_setpoint q_3_setpoint]);


vehicle_local_position=log.data.vehicle_local_position_0{:,:}; 
XYZ=vehicle_local_position(:,6:8);
V_XYZ=vehicle_local_position(:,12:14);

vehicle_local_position_setpoint=log.data.vehicle_local_position_setpoint_0{:,:}; 
XYZ_setpoint=vehicle_local_position_setpoint(:,2:4);
V_XYZ_setpoint=vehicle_local_position_setpoint(:,7:9);

actuator_controls=log.data.actuator_controls_0_0{:,:};   
Roll_control=actuator_controls(:,3);
Pitch_control=actuator_controls(:,4);
Yaw_control=actuator_controls(:,5);


actuator_outputs=log.data.actuator_outputs_0{:,:};   
cs1=actuator_outputs(:,5);
cs2=actuator_outputs(:,6);
cs3=actuator_outputs(:,7);
cs4=actuator_outputs(:,8);


indi_feedback_input=log.data.indi_feedback_input_0{:,:};     


[rate_N,~]=size(vehicle_rates_setpoint(:,1));
rate_delta_t=zeros(rate_N-1,1);
for i=1:rate_N-1
rate_delta_t(i)=(vehicle_rates_setpoint(i+1,1))*1e-6-(vehicle_rates_setpoint(i,1))*1e-6;
end

[attitude_N,~]=size(vehicle_attitude_setpoint(:,1));
attitude_delta_t=zeros(attitude_N-1,1);
for i=1:attitude_N-1
attitude_delta_t(i)=(vehicle_attitude_setpoint(i+1,1))*1e-6-(vehicle_attitude_setpoint(i,1))*1e-6;
end

[pose_N,~]=size(vehicle_local_position_setpoint(:,1));
pose_delta_t=zeros(pose_N-1,1);
for i=1:pose_N-1
pose_delta_t(i)=(vehicle_local_position_setpoint(i+1,1))*1e-6-(vehicle_local_position_setpoint(i,1))*1e-6;
end

[actuator_N,~]=size(actuator_controls(:,1));
actuator_delta_t=zeros(actuator_N-1,1);
for i=1:actuator_N-1
actuator_delta_t(i)=(actuator_controls(i+1,1))*1e-6-(actuator_controls(i,1))*1e-6;
end

[cs_N,~]=size(actuator_outputs(:,1));
cs_delta_t=zeros(cs_N-1,1);
for i=1:cs_N-1
cs_delta_t(i)=(actuator_outputs(i+1,1))*1e-6-(actuator_outputs(i,1))*1e-6;
end
%% 
figure,
subplot(311)
plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,2)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,3)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('p(deg/s)')
title('Angular velocity');
legend('Setpoint','Response');
%% 
subplot(312)
plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,3)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,4)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('q(deg/s)')
% title('Pitch angular velocity');
legend('Setpoint','Response');
%% 
subplot(313)
plot((vehicle_rates_setpoint(:,1))*1e-6, vehicle_rates_setpoint(:,4)*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_angular_velocity(:,1))*1e-6, vehicle_angular_velocity(:,5)*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('r(deg/s)')
% title('Yaw angular velocity');
legend('Setpoint','Response');
%% 




%% and maybe more figure, all in the variable "log.data"
figure,
subplot(311)
plot((vehicle_attitude_setpoint(:,1))*1e-6, Roll_setpoint*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_attitude(:,1))*1e-6, Roll*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Roll(deg)')
title('Euler angular');
legend('Setpoint','Response');
%% and maybe more figure, all in the variable "log.data"
subplot(312)
plot((vehicle_attitude_setpoint(:,1))*1e-6, Pitch_setpoint*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_attitude(:,1))*1e-6, Pitch*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Pitch(deg)')
legend('Setpoint','Response');
%% and maybe more figure, all in the variable "log.data"
subplot(313)
plot((vehicle_attitude_setpoint(:,1))*1e-6, Yaw_setpoint*r2d,'k-','LineWidth',1);hold on;
plot((vehicle_attitude(:,1))*1e-6, Yaw*r2d,'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Yaw(deg)')
legend('Setpoint','Response');
%% 






%% and maybe more figure, all in the variable "log.data"
figure,
subplot(311)
plot((vehicle_local_position_setpoint(:,1))*1e-6, V_XYZ_setpoint(:,1),'k-','LineWidth',1);hold on;
plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,1),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
title('Velocity');
xlabel({'Time(s)'});
ylabel('V_X(m/s)')
legend('Setpoint','Response');
%% and maybe more figure, all in the variable "log.data"
subplot(312)
plot((vehicle_local_position_setpoint(:,1))*1e-6, V_XYZ_setpoint(:,2),'k-','LineWidth',1);hold on;
plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,2),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('V_Y(m/s)')
legend('Setpoint','Response');
%% and maybe more figure, all in the variable "log.data"
subplot(313)
plot((vehicle_local_position_setpoint(:,1))*1e-6, V_XYZ_setpoint(:,3),'k-','LineWidth',1);hold on;
plot((vehicle_local_position(:,1))*1e-6, V_XYZ(:,3),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('V_Z(m/s)')
legend('Setpoint','Response');
%% 





%% and maybe more figure, all in the variable "log.data"
figure,
subplot(311)
plot((vehicle_local_position_setpoint(:,1))*1e-6, XYZ_setpoint(:,1),'k-','LineWidth',1);hold on;
plot((vehicle_local_position(:,1))*1e-6, XYZ(:,1),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
title('Position');
xlabel({'Time(s)'});
ylabel('X(m)')
legend('Setpoint','Response');
%% and maybe more figure, all in the variable "log.data"
subplot(312)
plot((vehicle_local_position_setpoint(:,1))*1e-6, XYZ_setpoint(:,2),'k-','LineWidth',1);hold on;
plot((vehicle_local_position(:,1))*1e-6, XYZ(:,2),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Y(m)')
legend('Setpoint','Response');
%% and maybe more figure, all in the variable "log.data"
subplot(313)
plot((vehicle_local_position_setpoint(:,1))*1e-6, XYZ_setpoint(:,3),'k-','LineWidth',1);hold on;
plot((vehicle_local_position(:,1))*1e-6, XYZ(:,3),'--','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
% axis([-inf inf -100 100]);
xlabel({'Time(s)'});
ylabel('Z(m)')
legend('Setpoint','Response');
%% 





%% 
figure,
subplot(411)
plot((actuator_controls(:,1))*1e-6, Roll_control(:,1),'r-','LineWidth',1);hold on;
plot((actuator_controls(:,1))*1e-6, Pitch_control(:,1),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
plot((actuator_controls(:,1))*1e-6, Yaw_control(:,1),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
title('Actuator controls');
xlabel({'Time(s)'});
ylabel('Control')
legend('Roll','Pitch',  'Yaw');
subplot(412)
plot((actuator_controls(:,1))*1e-6, Roll_control(:,1),'r-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Roll control')
subplot(413)
plot((actuator_controls(:,1))*1e-6, Pitch_control(:,1),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Pitch control')
subplot(414)
plot((actuator_controls(:,1))*1e-6, Yaw_control(:,1),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Yaw control')
%% 





%% 
figure,
subplot(411)
plot((indi_feedback_input(:,1))*1e-6, indi_feedback_input(:,3),'r-','LineWidth',1);hold on;
plot((indi_feedback_input(:,1))*1e-6, indi_feedback_input(:,4),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
plot((indi_feedback_input(:,1))*1e-6, indi_feedback_input(:,5),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
title('indi feedback input');
xlabel({'Time(s)'});
ylabel('Control')
legend('Roll','Pitch', 'Yaw');
subplot(412)
plot((indi_feedback_input(:,1))*1e-6, indi_feedback_input(:,3),'r-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Roll control')
subplot(413)
plot((indi_feedback_input(:,1))*1e-6, indi_feedback_input(:,4),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Pitch control')
subplot(414)
plot((indi_feedback_input(:,1))*1e-6, indi_feedback_input(:,5),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Yaw control')
%% 


%% 
figure,
subplot(411)
plot((actuator_controls(:,1))*1e-6, Roll_control(:,1)-indi_feedback_input(:,3),'r-','LineWidth',1);hold on;
plot((actuator_controls(:,1))*1e-6, Pitch_control(:,1)-indi_feedback_input(:,4),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
plot((actuator_controls(:,1))*1e-6, Yaw_control(:,1)-indi_feedback_input(:,5),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
title('dynamic controls');
xlabel({'Time(s)'});
ylabel('Control')
legend('Roll','Pitch',  'Yaw');
subplot(412)
plot((actuator_controls(:,1))*1e-6, Roll_control(:,1)-indi_feedback_input(:,3),'r-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Roll control')
subplot(413)
plot((actuator_controls(:,1))*1e-6, Pitch_control(:,1)-indi_feedback_input(:,4),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Pitch control')
subplot(414)
plot((actuator_controls(:,1))*1e-6, Yaw_control(:,1)-indi_feedback_input(:,5),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf -0.5 0.5]);
xlabel({'Time(s)'});
ylabel('Yaw control')
%% 


%% 
figure,
subplot(511)
plot((actuator_outputs(:,1))*1e-6, cs1(:,1),'r-','LineWidth',1);hold on;
plot((actuator_outputs(:,1))*1e-6, cs2(:,1),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
plot((actuator_outputs(:,1))*1e-6, cs3(:,1),'b-','LineWidth',1);hold on;
plot((actuator_outputs(:,1))*1e-6, cs4(:,1),'g-','LineWidth',1);hold on;
grid on;
axis([-inf inf 800 2000]);
title('control surface');
xlabel({'Time(s)'});
ylabel('Deflection(pwm)')
legend('cs1','cs2','cs3','cs4');
subplot(512)
plot((actuator_outputs(:,1))*1e-6, cs1(:,1),'r-','LineWidth',1);hold on;
grid on;
axis([-inf inf 800 2000]);
xlabel({'Time(s)'});
ylabel('cs1')
subplot(513)
plot((actuator_outputs(:,1))*1e-6, cs2(:,1),'k-','LineWidth',1,'color',[0.6,0.2,0]);hold on;
grid on;
axis([-inf inf 800 2000]);
xlabel({'Time(s)'});
ylabel('cs2')
subplot(514)
plot((actuator_outputs(:,1))*1e-6, cs3(:,1),'b-','LineWidth',1);hold on;
grid on;
axis([-inf inf 800 2000]);
xlabel({'Time(s)'});
ylabel('cs3')
subplot(515)
plot((actuator_outputs(:,1))*1e-6, cs4(:,1),'g-','LineWidth',1);hold on;
grid on;
axis([-inf inf 800 2000]);
xlabel({'Time(s)'});
ylabel('cs4')
%% 

% 
% %% 
% figure,
% subplot(311)
% plot(1:rate_N-1, rate_delta_t,'k-','LineWidth',1);hold on;
% subplot(312)
% plot(1:attitude_N-1, attitude_delta_t,'k-','LineWidth',1);hold on;
% subplot(313)
% plot(1:pose_N-1, pose_delta_t,'k-','LineWidth',1);hold on;
