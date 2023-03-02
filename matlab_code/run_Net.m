function run_Net(N,net)

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
% [networks2.(fields{net}), sim_data.(fields{net})] = simulate_ms_sys(smallNetworks.(fields{net}),U);
% [networks2.(fields{net}), sim_data.(fields{net})] = simulate_hysteresis_sys(smallNetworks.(fields{net}),U);
[networks2.(fields{net}), sim_data.(fields{net})] = simulate_maxwell_sys(smallNetworks.(fields{net}), U);

save(['files/smallNetworksMemory_maxwell/networks2_net_' num2str(net) '_N_' num2str(N)],'networks2');
save(['files/smallNetworksMemory_maxwell/sim_data_net_' num2str(net) '_N_' num2str(N)],'sim_data');


end