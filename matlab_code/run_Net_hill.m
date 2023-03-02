function run_Net_hill(N,net)

% Load in the desried networks
load('files/smallNetworks4.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data = learningTask(N,250000); % I/O supervised learning task

wash_out = 60000;
% wash_out = 40000;
start = 80000-wash_out;
%start = 20000;
len = 200000+start;
%len = 80000+start;
len_test = 15000;

Data.u = (Data.u-min(Data.u))./(max(Data.u)-min(Data.u));
U = Data.u(start:len,1);
Y = Data.y(start:len,1);

% prepare testing data
U_test = Data.u(len+1:len+len_test,1);  
Y_test = Data.y(len+1:len+len_test,1); 


fields = fieldnames(smallNetworks);
% [networks2.(fields{net}), sim_data.(fields{net})] = simulate_ms_sys(smallNetworks.(fields{net}),U);
[networks2.(fields{net}), sim_data.(fields{net})] = simulate_hill_sys(smallNetworks.(fields{net}),U);
% [networks2.(fields{net}), sim_data.(fields{net})] = simulate_maxwell_sys(smallNetworks.(fields{net}), U);

save(['files/smallNets_hill_activation/networks2_hill_net_' num2str(net) '_N_' num2str(N)],'networks2');
save(['files/smallNets_hill_activation/sim_data_hill_net_' num2str(net) '_N_' num2str(N)],'sim_data');


end
