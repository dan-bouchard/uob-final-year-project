function run_Net_hill_complexity(N,net)

% Load in the desried networks
load('files/smallNetworks4.mat');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load(['files/complexity_data/Complexity' num2str(N)]);

wash_out = 60000;
start = 80000-wash_out;
len = 200000+start;
len_test = 15000;


dat.u = (dat.u-min(dat.u))./(max(dat.u)-min(dat.u));
U = dat.u(start:len,1);  
Y = dat.y(start:len,1);

% prepare testing data
U_test = dat.u(len+1:len+len_test,1);  
Y_test = dat.y(len+1:len+len_test,1); 


fields = fieldnames(smallNetworks);
% [networks2.(fields{net}), sim_data.(fields{net})] = simulate_ms_sys(smallNetworks.(fields{net}),U);
[networks2.(fields{net}), sim_data.(fields{net})] = simulate_hill_sys(smallNetworks.(fields{net}),U);
% [networks2.(fields{net}), sim_data.(fields{net})] = simulate_maxwell_sys(smallNetworks.(fields{net}), U);

save(['files/smallNets_complexity_hill/networks2_hill_net_' num2str(net) '_N_' num2str(N)],'networks2');
save(['files/smallNets_complexity_hill/sim_data_hill_net_' num2str(net) '_N_' num2str(N)],'sim_data');


end
