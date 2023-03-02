function [d] = simulate_one_network_testA(N)
% Run the whole simulation for only one network

load('files/smallNetworks2.mat');
net = smallNetworks.net7;
net.init_data.save_sim_data = 1;


Data = learningTask(N,250000); % I/O supervised learning task

wash_out = 60000;
% start = 80000-wash_out;
start = 20000;
len = 200000+start;

len_test = 15000;

U = Data.u(start:len,1);  
Y = Data.y(start:len,1);

% prepare testing data
U_test = Data.u(len+1:len+len_test,1);  
Y_test = Data.y(len+1:len+len_test,1); 



% [networks2, sim_data] = simulate_hysteresis_sys(net, U);
[networks2, sim_data] = simulate_ms_sys(net,U);
% [networks2, sim_data] = simulate_maxwell_sys(net, U);
% [networks2, sim_data] = simulate_kelvin_a_sys(net, U);

% save('sim_data');

X = sim_data.D(wash_out:end,:);

Yw = Y(wash_out:end,:);

W_out=X\Yw; % X*W_out=Yw e.g. (1000x77)*(77x1)=(1000x1) 

net_test = networks2; %output network
net_test.W_out = W_out; % output weights for each spring

% [net_test_out,sim_data_test] = simulate_hysteresis_sys(net_test, U_test);
[net_test_out,sim_data_test] = simulate_ms_sys(net_test,U_test);
% [net_test_out,sim_data_test] = simulate_maxwell_sys(net_test,U_test);
% [net_test_out,sim_data_test] = simulate_kelvin_a_sys(net_test,U_test);
mse = mean_squared_error(Y_test,sim_data_test.O);

% % plot results
% figure;plot(Y_test,'r','LineWidth',1);
% hold on;plot(sim_data_test.O,'--','LineWidth',1);
% f1=gcf;a1=gca;
% set(a1,'FontSize',14);
% xlabel('timestep [ ]');
% ylabel('[ ]');
% title(['Performance comparison - Net 1 (MSE: ', num2str(mse), ')']);
% legend('target output','system output')

disp(['MSE: ',num2str(mse)]);
d = sim_data;

end