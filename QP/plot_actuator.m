% in=simout1.Data
% out=simout.Data;
% power_brushlessDCmotor
% clc;
clear all;
close all;
addpath D:\17master_mengchaoheng_data\qcat1_2_1\QCAT\qcat
% function callqpact
global B plim rlim umin umax vmin vmax a_c a_t com_t_wls com_t_inv com_t_dca dca_W2 com_c p_limits S dca_W1 select_test hf lf Am%Del
% Del=diag([20*pi/180 20*pi/180 20*pi/180 20*pi/180 20*pi/180 20*pi/180 20*pi/180 20*pi/180 200],0);%
% 转速差作为新操纵面%论文数据
% p_limits=[20*pi/180 200];%[rad  rad/s]
% v_limits=[400*pi/180 800];% [rad/s rad/s^2]
% 转速差作为新操纵面%实测数据
% p_limits=[20*pi/180 200];%[rad  rad/s]
% v_limits=[400*pi/180 466.6667];% [rad/s rad/s^2]
% 拉力差作为新操纵面%论文数据
% p_limits=[20*pi/180 4.89];%[rad  rad/s]
% v_limits=[400*pi/180 19.56];% [rad/s rad/s^2]
% 拉力差作为新操纵面%实测数据
p_limits=[20*pi/180 4.2513];%[rad  rad/s]% 2*k_TS*speed*200=4.2513 % 200为用于提供滚转力矩的转速差+-200
v_limits=[400*pi/180 9.9197];% [rad/s rad/s^2]% 实测转速最大变化率466.6667rad/s^2，2*k_TS*speed*466.6667=9.9197

%===============================================================
k_TS=9.9796018325697625989171178675552e-6;
speed=1065;
kc=3.157;
l_1=0.17078793-0.09;% roll,pitch
l_2=0.175;% T
l_3=0.06647954;% yaw1 yaw3
l_4=0.06647954+0.175;% yaw4
l_5=0.175-0.06647954;% yaw2
I_x=0.054593;
I_y=0.017045;
I_z=0.049226;
I=[I_x 0 0;0 I_y 0;0 0 I_z];
% =========================9 sureface=================================
L=[-l_1 0 l_1 0 -l_1 0 l_1 0 l_2;
    0 -l_1 0 l_1 0 -l_1 0 l_1 0;
    l_3 -l_5 l_3 l_4 l_3 l_4 l_3 -l_5 0];
% F=diag([kc kc kc kc kc kc kc kc 2*k_TS*speed],0);%% 转速差作为新操纵面
F=diag([kc kc kc kc kc kc kc kc 1],0);% 拉力差作为新操纵面
B=I\L*F;
[m,k] = size(B);
a_c=35;
a_t=4;
%0.5%6%Am=6
%0.5%5%Am=8
% lf=0.5;hf=3;Am=8;
lf=0.5;hf=6;Am=6;
W1_c=p_limits(1);
W1_t=p_limits(2);
W2_c=2*p_limits(1)/v_limits(1);
W2_t=2*p_limits(2)/v_limits(2);
% dca_W1=diag([W1_c W1_c W1_c W1_c W1_c W1_c W1_c W1_c W1_t]);% eye(9);%
dca_W2=1*diag([W2_c W2_c W2_c W2_c W2_c W2_c W2_c W2_c W2_t]);
dca_W1=1*eye(9);% 
% dca_W2=zeros(9,9);
% dca_W2=diag([1 1 1 1 1 1 1 1 100]);
% 频率信号
select_test=0;
com_c=1/(1-exp(-90*0.01));% 90
com_t_wls=1/(1-exp(-10*0.01));%11
com_t_dca=1/(1-exp(-10*0.01));
com_t_inv=1/(1-exp(-10*0.01));

% com_c=1;
% com_t_wls=1;
% com_t_dca=1;
% com_t_inv=1;
% load('Actuator.mat');
S=zeros(k,m);
P_1=pinv([B(:,2) B(:,4) B(:,6) B(:,8) B(:,9)])
% P_1=pinv([B(:,4) B(:,6) B(:,9)]);
% P_1=pinv(B(:,1:8));
% S=pinv(B)
S([2 4 6 8 9],:)=P_1;
% S(1:8,:)=P_1;
% S([4 6 9],:)=P_1;
%========================================
% S=zeros(k,m);

% Ws=diag([1 1 1 1 0.01]);
% Ws=diag([1 1 0.01]);
% B1=[B(:,2) B(:,4) B(:,6) B(:,8) B(:,9)];
% B1=[B(:,4) B(:,6) B(:,9)];
% S([2 4 6 8 9],:)=Ws\B1'/( B1*(Ws\B1') );
% (Ws\B1')/( B1*(Ws\B1') )
% B1'/( B1*B1' )
% pinv(B1)
% Ws=diag([1/W1_c 1/W1_c 1/W1_c 1/W1_c 1/W1_c 1/W1_c 1/W1_c 1/W1_c 1/W1_t]);
% S=Ws\B'/( B*(Ws\B') )
% pinv(B)
%========================================
umin=[[1;1;1;1;1;1;1;1]*(-p_limits(1));-p_limits(2)];
umax=[[1;1;1;1;1;1;1;1]*p_limits(1);p_limits(2)];
vmin=[[1;1;1;1;1;1;1;1]*(-v_limits(1));-v_limits(2)];
vmax=[[1;1;1;1;1;1;1;1]*v_limits(1);v_limits(2)];
% vmin=Del\[[1;1;1;1;1;1;1;1]*(-v_limits(1));-v_limits(2)];
% vmax=Del\[[1;1;1;1;1;1;1;1]*v_limits(1);v_limits(2)];
plim=[umin umax];
rlim=[vmin vmax];% [vmin vmax];
sim('twin');%一、0.1-8HZ，幅值5.二、阶跃6
in_x=input.Data(:,1);
in_y=input.Data(:,2);
in_z=input.Data(:,3);
out_wls_x_off=output_wls_off.Data(:,1);
out_wls_y_off=output_wls_off.Data(:,2);
out_wls_z_off=output_wls_off.Data(:,3);

out_wls_x=output_wls.Data(:,1);
out_wls_y=output_wls.Data(:,2);
out_wls_z=output_wls.Data(:,3);

out_dca_x=output_dca.Data(:,1);
out_dca_y=output_dca.Data(:,2);
out_dca_z=output_dca.Data(:,3);

dca_surfase_in=dca_c_in.Data;
dca_surfase_out=dca_c_out.Data;
dca_T_in=dca_in.Data;
dca_T_out=dca_out.Data;

wls_surfase_in=wls_c_in.Data;
wls_surfase_out=wls_c_out.Data;
wls_T_in=wls_in.Data;
wls_T_out=wls_out.Data;

wls_surfase_in_off=wls_c_in_off.Data;
wls_surfase_out_off=wls_c_out_off.Data;
wls_T_in_off=wls_in_off.Data;
wls_T_out_off=wls_out_off.Data;
% inv_surfase_in=inv_c_in.Data;
% inv_surfase_out=inv_c_out.Data;
% inv_T_in=inv_in.Data;
% inv_T_out=inv_out.Data;

time=0:0.01:3;
tt=1:30:301;
% figure,%
% subplot(3,1,1)
% plot(time,in_x,'r-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_x,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_dca_x,'Color','g','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_x,'Color','k','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % axis([0 2 -20 20]);
% title('伪控制指令响应曲线','FontSize',8);
% xlabel('\itt \rm(s)','FontSize',8);ylabel('\it\Gamma_p \rm(N*m)','FontSize',8)
% h1=legend('cmd','dca','wls');% ,'Location','EastOutside';,'\tau_{i}'
% set(h1,'FontName','Times New Roman','FontSize',8,'NumColumns',4,'location','southwest');
% 
% subplot(3,1,2)
% % figure,%
% plot(time,in_y,'r-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_y,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_dca_y,'Color','g','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_y,'Color','k','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % axis([0 2 -20 20]);
% xlabel('\itt \rm(s)','FontSize',8);ylabel('\it\Gamma_q \rm(N*m)','FontSize',8)
% 
% subplot(3,1,3)
% % figure,%
% plot(time,in_z,'r-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_z,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_dca_z,'Color','g','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_z,'Color','k','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % axis([0 2 -20 20]);
% xlabel('\itt \rm(s)','FontSize',8);ylabel('\it\Gamma_r \rm(N*m)','FontSize',8)

figure,%
subplot(3,1,1)
% plot(time,dca_surfase_in,'b-','Marker','o','MarkerIndices',tt);hold on;
plot(time,dca_surfase_out,'b-','Marker','none','MarkerIndices',tt);hold on;
% plot(time,inv_surfase_in/(p_limits(1)),'g-','Marker','o','MarkerIndices',tt);hold on;
% plot(time,inv_surfase_out/(p_limits(1)),'g--','Marker','*','MarkerIndices',tt);hold on;
% plot(time,wls_surfase_in,'r-','Marker','s','MarkerIndices',tt);hold on;
plot(time,wls_surfase_out,'r--','Marker','none','MarkerIndices',tt);hold on;grid on;
% axis([0 3 -0.52 0.35]);
xlabel('Time (s)');
ylabel('Deflection (rad)')
title('Control effector \delta_3');
% h2=legend('cmd_3^{d}','u_3^{d}','cmd_3^{s}','u_3^{s}');% ,'cmd_1^{i}','u_1^{i}'
h2=legend('dca','sca');
% h2=legend('\delta^d','\delta^s');
set(h2,'NumColumns',1,'location','southwest');%northwest
set(gca, 'FontSize', 10,'FontName','Times New Roman')
% figure,%
subplot(3,1,2)
% plot(time,dca_T_in,'b-','Marker','o','MarkerIndices',tt);hold on;
plot(time,dca_T_out,'b-','Marker','none','MarkerIndices',tt);hold on;
% plot(time,inv_T_in/(p_limits(2)),'g-','Marker','o','MarkerIndices',tt);hold on;
% plot(time,inv_T_out/(p_limits(2)),'g--','Marker','*','MarkerIndices',tt);hold on;
% plot(time,wls_T_in,'r','Marker','s','MarkerIndices',tt);hold on;
plot(time,wls_T_out,'r--','Marker','none','MarkerIndices',tt);hold on;grid on;
% axis([0 3 -100 85]);
xlabel('Time (s)');
ylabel('Speed difference (rad/s)');
title('Control effector \delta_9');
% set(gca,'FontName','Times New Roman');
% h3=legend('cmd_9^{d}','u_9^{d}','cmd_9^{s}','u_9^{s}');% ,'cmd_9^{i}','u_9^{i}'
h3=legend('dca','sca');
% h3=legend('\delta^d','\delta^s');
set(h3,'NumColumns',1,'location','southwest');
set(gca, 'FontSize', 10,'FontName','Times New Roman')%
% figure,%
subplot(3,1,3)
plot(time,in_x,'k-','Marker','none','MarkerIndices',tt);hold on;
% plot(time,out_inv_x,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
plot(time,out_dca_x,'Color','b','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
plot(time,out_wls_x,'Color','r','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% axis([0 3 -10.5 6.6]);
title('Virtual input');
xlabel('Time (s)','FontName','Times New Roman');
ylabel('$\dot p\ \rm(rad/s^2)$','interpreter','latex')
h1=legend('cmd','dca','sca');% ,'Location','EastOutside';,'\tau_{i}'
set(h1,'NumColumns',1,'location','southwest');
set(gca, 'FontSize', 10,'FontName','Times New Roman')
%------------------------------------------------------------------------------------------------------------------------
% 阶跃信号
select_test=1;
com_c=1/(1-exp(-100*0.01));
com_t_wls=1/(1-exp(-11*0.01));
com_t_dca=1/(1-exp(-11*0.01));
sim('twin');%一、0.1-8HZ，幅值5.二、阶跃6
in_x=input.Data(:,1);
in_y=input.Data(:,2);
in_z=input.Data(:,3);
% out_inv_x=output_inv.Data(:,1);
% out_inv_y=output_inv.Data(:,2);
% out_inv_z=output_inv.Data(:,3);


out_wls_x=output_wls.Data(:,1);
out_wls_y=output_wls.Data(:,2);
out_wls_z=output_wls.Data(:,3);


out_dca_x=output_dca.Data(:,1);
out_dca_y=output_dca.Data(:,2);
out_dca_z=output_dca.Data(:,3);

dca_surfase_in=dca_c_in.Data;
dca_surfase_out=dca_c_out.Data;
dca_T_in=dca_in.Data;
dca_T_out=dca_out.Data;

wls_surfase_in=wls_c_in.Data;
wls_surfase_out=wls_c_out.Data;
wls_T_in=wls_in.Data;
wls_T_out=wls_out.Data;

% inv_surfase_in=inv_c_in.Data;
% inv_surfase_out=inv_c_out.Data;
% inv_T_in=inv_in.Data;
% inv_T_out=inv_out.Data;
tt=1:10:201;
% figure,%
% subplot(3,1,1)
% plot(time,in_x,'r-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_x,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_dca_x,'Color','g','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_x,'Color','k','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % axis([0 2 -20 20]);
% title('伪控制指令响应曲线','FontSize',8);
% xlabel('\itt \rm(s)','FontSize',8,'FontName','Times New Roman');ylabel('\it\Gamma_p \rm(N*m)','FontSize',8,'FontName','Times New Roman')
% h1=legend('\tau_{c}','\tau_{i}','\tau_{d}','\tau_{w}');% ,'Location','EastOutside';
% set(h1,'FontName','Times New Roman','FontSize',8,'NumColumns',3,'location','East');
% 
% subplot(3,1,2)
% % figure,%
% plot(time,in_y,'r-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_y,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_dca_y,'Color','g','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_y,'Color','k','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % axis([0 2 -20 20]);
% xlabel('\itt \rm(s)','FontSize',8,'FontName','Times New Roman');ylabel('\it\Gamma_q \rm(N*m)','FontSize',8,'FontName','Times New Roman')
% 
% subplot(3,1,3)
% % figure,%
% plot(time,in_z,'r-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_z,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_dca_z,'Color','g','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_z,'Color','k','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% % axis([0 2 -20 20]);
% xlabel('\itt \rm(s)','FontSize',8,'FontName','Times New Roman');ylabel('\it\Gamma_r \rm(N*m)','FontSize',8,'FontName','Times New Roman')

figure,%
subplot(3,1,1)
% plot(time,dca_surfase_in,'b-','Marker','o','MarkerIndices',tt);hold on;
plot(time,dca_surfase_out,'b-','Marker','none','MarkerIndices',tt);hold on;
% plot(time,inv_surfase_in/(p_limits(1)),'g-','Marker','o','MarkerIndices',tt);hold on;
% plot(time,inv_surfase_out/(p_limits(1)),'g--','Marker','*','MarkerIndices',tt);hold on;
% plot(time,wls_surfase_in,'r-','Marker','s','MarkerIndices',tt);hold on;
plot(time,wls_surfase_out,'r--','Marker','none','MarkerIndices',tt);hold on;grid on;
axis([0 0.5 -0.15 0.45]);
xlabel('Time (s)','FontSize',8,'FontName','Times New Roman');
ylabel('Deflection (rad)','FontSize',8,'FontName','Times New Roman')
title('Control effector \delta_3','FontSize',8);
% h2=legend('cmd_3^{d}','u_3^{d}','cmd_3^{s}','u_3^{s}');% ,'cmd_1^{i}','u_1^{i}'
h2=legend('dca','sca');
% h2=legend('\delta^d','\delta^s');
set(h2,'FontName','Times New Roman','FontSize',8,'NumColumns',1,'location','northwest');%northwest

subplot(3,1,2)
% figure,%
% plot(time,dca_T_in,'b-','Marker','o','MarkerIndices',tt);hold on;
plot(time,dca_T_out,'b-','Marker','none','MarkerIndices',tt);hold on;
% plot(time,inv_T_in/(p_limits(2)),'g-','Marker','o','MarkerIndices',tt);hold on;
% plot(time,inv_T_out/(p_limits(2)),'g--','Marker','*','MarkerIndices',tt);hold on;
% plot(time,wls_T_in,'r','Marker','s','MarkerIndices',tt);hold on;
plot(time,wls_T_out,'r--','Marker','none','MarkerIndices',tt);hold on;grid on;
axis([0 0.5 -10 90]);
xlabel('Time (s)','FontSize',8,'FontName','Times New Roman');
ylabel('Speed difference (rad/s)','FontSize',8,'FontName','Times New Roman');
title('Control effector \delta_9','FontSize',8);
% h3=legend('cmd_9^{d}','u_9^{d}','cmd_9^{s}','u_9^{s}');% ,'cmd_9^{i}','u_9^{i}'
h3=legend('dca','sca');
% h3=legend('\delta^d','\delta^s');
set(h3,'FontName','Times New Roman','FontSize',8,'NumColumns',1,'location','northwest');

subplot(3,1,3)
plot(time,in_x,'k-','Marker','none','MarkerIndices',tt);hold on;
% plot(time,out_inv_x,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
plot(time,out_dca_x,'Color','b','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
plot(time,out_wls_x,'Color','r','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
axis([0 0.5 -1 7]);
title('Virtual input','FontSize',8);
xlabel('Time (s)','FontSize',8,'FontName','Times New Roman');
ylabel('$\dot p\ \rm(rad/s^2)$','interpreter','latex','FontSize',8)
% h1=legend('cmd','dca','sca');% ,'Location','EastOutside';,'\tau_{i}'
h1=legend('cmd','dca','sca');
set(h1,'FontName','Times New Roman','FontSize',8,'NumColumns',1,'location','northwest');% southwest
% --------------执行器补偿
% 
% select_test=0;
% hf=10;
% com_c=1/(1-exp(-80*0.01));
% com_t_wls=1/(1-exp(-10*0.01));
% % com_c=1;
% % com_t_wls=1;
% sim('twin');%一、0.1-8HZ，幅值5.二、阶跃6
% in_x=input.Data(:,1);
% in_y=input.Data(:,2);
% in_z=input.Data(:,3);
% out_wls_x_off=output_wls_off.Data(:,1);
% out_wls_y_off=output_wls_off.Data(:,2);
% out_wls_z_off=output_wls_off.Data(:,3);
% wls_surfase_in_off=wls_c_in_off.Data;
% wls_surfase_out_off=wls_c_out_off.Data;
% wls_T_in_off=wls_in_off.Data;
% wls_T_out_off=wls_out_off.Data;
% 
% out_wls_x=output_wls.Data(:,1);
% out_wls_y=output_wls.Data(:,2);
% out_wls_z=output_wls.Data(:,3);
% wls_surfase_in=wls_c_in.Data;
% wls_surfase_out=wls_c_out.Data;
% wls_T_in=wls_in.Data;
% wls_T_out=wls_out.Data;
% 
% 
% 
% figure,%
% subplot(3,1,1)
% plot(time,wls_surfase_in_off,'k-','Marker','o','MarkerIndices',tt);hold on;
% plot(time,wls_surfase_out_off,'b-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,inv_surfase_in/(p_limits(1)),'g-','Marker','o','MarkerIndices',tt);hold on;
% % plot(time,inv_surfase_out/(p_limits(1)),'g--','Marker','*','MarkerIndices',tt);hold on;
% plot(time,wls_surfase_in,'k-','Marker','s','MarkerIndices',tt);hold on;
% plot(time,wls_surfase_out,'r--','Marker','none','MarkerIndices',tt);hold on;grid on;
% % axis([0 3 -1.3 1]);
% xlabel('\itt \rm(s)','FontSize',8,'FontName','Times New Roman');ylabel('\itu_3 \rm(rad)','FontSize',8,'FontName','Times New Roman')
% % h2=legend('cmd_3^{d}','u_3^{d}','cmd_3^{s}','u_3^{s}');% ,'cmd_1^{i}','u_1^{i}'
% h2=legend('cmd','dynamic','compensate');;
% set(h2,'FontName','Times New Roman','FontSize',8,'NumColumns',1,'location','best');%northwest
% 
% % figure,%
% subplot(3,1,2)
% plot(time,wls_T_in_off,'k-','Marker','o','MarkerIndices',tt);hold on;
% plot(time,wls_T_out_off,'b-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,inv_T_in/(p_limits(2)),'g-','Marker','o','MarkerIndices',tt);hold on;
% % plot(time,inv_T_out/(p_limits(2)),'g--','Marker','*','MarkerIndices',tt);hold on;
% plot(time,wls_T_in,'k','Marker','s','MarkerIndices',tt);hold on;
% plot(time,wls_T_out,'r--','Marker','none','MarkerIndices',tt);hold on;grid on;
% % axis([0 3 -0.55 0.45]);
% xlabel('\itt \rm(s)','FontSize',8,'FontName','Times New Roman');
% ylabel('\itu_9 \rm(rad/s)','FontSize',8,'FontName','Times New Roman');
% % set(gca,'FontName','Times New Roman');
% % h3=legend('cmd_9^{d}','u_9^{d}','cmd_9^{s}','u_9^{s}');% ,'cmd_9^{i}','u_9^{i}'
% % h3=legend('cmd','dynamic','compensate');;
% set(h3,'FontName','Times New Roman','FontSize',8,'NumColumns',3,'location','best');
% 
% % figure,%
% subplot(3,1,3)
% plot(time,in_x,'k-','Marker','none','MarkerIndices',tt);hold on;
% % plot(time,out_inv_x,'Color','b','LineStyle','-','Marker','^','MarkerIndices',tt);hold on;
% plot(time,out_wls_x_off,'Color','b','LineStyle','-.','Marker','x','MarkerIndices',tt);hold on;grid on;
% plot(time,out_wls_x,'Color','r','LineStyle','--','Marker','o','MarkerIndices',tt);hold on;
% axis([0 3 -10.5 6.6]);
% % title('伪控制指令响应曲线','FontSize',8);
% xlabel('\itt \rm(s)','FontSize',8,'FontName','Times New Roman');
% ylabel('$\dot p\ \rm(rad/s^2)$','interpreter','latex','FontSize',8)
% h1=legend('cmd','dynamic','compensate');% ,'Location','EastOutside';,'\tau_{i}'
% set(h1,'FontName','Times New Roman','FontSize',8,'NumColumns',3,'location','best');