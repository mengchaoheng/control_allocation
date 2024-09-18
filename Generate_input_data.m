clear all;
close all;
addpath(genpath(pwd))

% run 'plot_fly_log_states.m' for command_px4 and u_px4
load 'flight.mat'; % command_px4 and u_px4
[len_command_px4,~]=size(command_px4);
M=50;
[X,Y,Z] = sphere(M);
unit_vector=zeros(3,(M+1)^2);
N=(M+1)^2;

for i=1:N% (N+1)^2  for  sphere %length(M_des(1:1000,1))%%length(X)
    unit_vector(:,i)=200*[X(i); Y(i); Z(i)];
end
v=[command_px4'  zeros(3,200) unit_vector];% zeros is just for test

% 指定保存路径及文件名
filename = 'input.mat';
% 调用 save 函数进行保存
save(filename, 'v',"len_command_px4",'u_px4','delta_t_s');

filename = 'input.csv';
writematrix(v',filename);