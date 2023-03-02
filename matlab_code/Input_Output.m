function Input_Output(N)

% Get the input output graph out

load(['files/smallNetworksMemory/networks2_N_' num2str(N)]);
load(['files/smallNetworksMemory/sim_data_N_' num2str(N)]);

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

X = sim_data.net3.D(wash_out:end,:);
Yw = Y(wash_out:end,:);

W_out=X\Yw

net_test = networks2.net3; %output network
net_test.W_out = W_out; % output weights for each spring

[net_test_out,sim_data_test] = simulate_ms_sys(net_test,U_test);

% plot results
figure;plot(Y_test,'r','LineWidth',1);
hold on;plot(sim_data_test.O,'--','LineWidth',1);
f1=gcf;a1=gca;
set(a1,'FontSize',14);
xlabel('timestep [ ]');
ylabel('[ ]');
title('Performance comparison - Net 3')
legend('target output','system output')

disp(['MSE: ',num2str(mean_squared_error(Y_test,sim_data_test.O))]);
end