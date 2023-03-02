function [] = calculateResultsOneNet(N, network_num)
% N is the amount of memory of the task defined
% e.g. N=2
% must run runOneNetwork first before calculating results



load('files/smallNetworks3.mat');
load(['files/smallNets3_Test/networks2_hys_net_' num2str(network_num) '_net_' num2str(network_num)]);
load(['files/TEST/sim_data_N_' num2str(N) '_net_' num2str(network_num)]);

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

% Error = zeros(10,1);
net_fields = fieldnames(smallNetworks);


if (strcmp(smallNetworks.(net_fields{network_num}).readout_type,'LENGTHS'))
    % sim_data.D lengths of every spring
    X = sim_data.D(wash_out:end,:);  % throw first 100 steps away
else
    % sim_data.Sx x-positions of each mass
    X = sim_data.Sx(wash_out:end,:);  % throw first 100 steps away
end

Yw = Y(wash_out:end,:);


W_out=X\Yw; % X*W_out=Yw e.g. (1000x77)*(77x1)=(1000x1) 


net_test = networks2; %output network
net_test.W_out = W_out; % output weights for each spring


[net_test_out,sim_data_test] = simulate_ms_sys(net_test,U_test);
Error = mean_squared_error(Y_test,sim_data_test.O);
disp(['MSE: ',num2str(Error)])


save(['files/smallNets3_Test/Error_N_' num2str(N) '_net_' num2str(network_num)],'Error');

figure;plot(Y_test,'r','LineWidth',1);
hold on;plot((mapstd(sim_data_test.O'))','--','LineWidth',1);
legend('target output','system output');


