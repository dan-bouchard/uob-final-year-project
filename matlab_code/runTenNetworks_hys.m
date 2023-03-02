function [] = runTenNetworks_hys(N)
% N is the amount of memory of the task defined
% e.g. N=2


% Load in the desried networks
load('files/smallNetworks3.mat');


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
for i=1:10
    disp([' j = ' ,num2str(i) , ' of  10']);
    [networks2.(fields{i}), sim_data.(fields{i})] = simulate_hysteresis_sys(smallNetworks.(fields{i}),U);
end

save(['files/smallNets3_Test/networks2_N_' num2str(N)],'networks2');
save(['files/smallNets3_Test/sim_data_N_' num2str(N)],'sim_data');
