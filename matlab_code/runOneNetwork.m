function [] = runOneNetwork(N, network_num)
% N is the amount of memory of the task defined
% e.g. N=2


% Load in the desried networks
load('files/smallNetworks.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = learningTask(N,250000); % I/O supervised learning task

wash_out = 60000;
start = 80000-wash_out; 
len = 200000+start;
len_test = 15000;

U = Data.u(start:len,1);  
Y = Data.y(start:len,1);

% prepare testing data
U_test = Data.u(len+1:len+len_test,1);  
Y_test = Data.y(len+1:len+len_test,1); 


fields = fieldnames(smallNetworks);
% network type is spring or hysteresis

[networks2, sim_data] = simulate_ms_sys(smallNetworks.(fields{network_num}),U);
% [networks2, sim_data] = simulate_hysteresis_sys(smallNetworks.(fields{1}),U);


save(['files/TEST/networks2_N_' num2str(N) '_net_' num2str(network_num)],'networks2');
save(['files/TEST/sim_data_N_' num2str(N) '_net_' num2str(network_num)],'sim_data');

clear all;
close all;