function [] = simulate_one_network(N)
% Run the whole simulation for only one network

load('files/smallNetworks4.mat');
net = smallNetworks.net3;
net.init_data.save_sim_data = 1;
% net.W.k1 = rand_in_range_exp([1 100],length(net.W.k1));
% net.W.k3 = rand_in_range_exp([0.1 10],length(net.W.k1));
% net.W.d1 = rand_in_range_exp([1 100],length(net.W.k1));
% net.W.d3 = rand_in_range_exp([0.1 10],length(net.W.k1));

% net.W.k1 = rand_in_range_exp([10 100],length(net.W.k1));
% net.W.k3 = rand_in_range_exp([4 60],length(net.W.k1));
% net.W.d1 = rand_in_range_exp([10 100],length(net.W.k1));
% net.W.d3 = rand_in_range_exp([4 60],length(net.W.k1));
% net.W.beta = 10*ones(size(net.W.k1));
% net.W.gamma = -40*ones(size(net.W.k1));
% net.W.beta = rand_in_range_exp([0.1 100],length(net.W.k1)).*(2*randi(2,length(net.W.k1),1)-3);
% net.W.gamma = rand_in_range_exp([0.1 100],length(net.W.k1)).*(2*randi(2,length(net.W.k1),1)-3);


Data = learningTask(N,260000); % I/O supervised learning task

wash_out = 60000;
%wash_out = 40000;
start = 80000-wash_out;
%start = 20000;
len = 200000+start;
%len = 80000+start;
len_test = 15000;

U = Data.u(start:len,1);  
Y = Data.y(start:len,1);

% prepare testing data
U_test = Data.u(len+1:len+len_test,1);  
Y_test = Data.y(len+1:len+len_test,1); 



%[networks2_hys, sim_data_hys] = simulate_hysteresis_sys(net, U);
[networks2, sim_data] = simulate_hill_sys(net,U);
% [networks2, sim_data] = simulate_maxwell_sys(net, U);
% [networks2, sim_data] = simulate_kelvin_a_sys(net, U);
% [networks2, sim_data] = simulate_kelvin_b_sys(net, U);

save(['files/save_sim_data/networks2_hill_N_' num2str(N)],'networks2');
save(['files/save_sim_data/sim_data_hill_N_' num2str(N)],'sim_data');

X = sim_data.D(wash_out:end,:);
%X_hys = sim_data_hys.D(wash_out:end,:);

Yw = Y(wash_out:end,:);


W_out=X\Yw; % X*W_out=Yw e.g. (1000x77)*(77x1)=(1000x1) 
%W_out_hys=X_hys\Yw;

net_test = networks2; %output network
%net_test_hys = networks2_hys; %output network
net_test.W_out = W_out; % output weights for each spring
%net_test_hys.W_out = W_out_hys; % output weights for each spring

%[net_test_out,sim_data_test_hys] = simulate_hysteresis_sys(net_test_hys, U_test);
[net_test_out,sim_data_test] = simulate_hill_sys(net_test,U_test);
% [net_test_out,sim_data_test] = simulate_maxwell_sys(net_test,U_test);
% [net_test_out,sim_data_test] = simulate_kelvin_a_sys(net_test,U_test);
% [net_test_out,sim_data_test] = simulate_kelvin_b_sys(net_test,U_test);

Error = mean_squared_error(Y_test,sim_data_test.O);
%mse_hys = mean_squared_error(Y_test,mapstd(sim_data_test_hys.O));

% plot results
% figure;plot(Y_test,'r','LineWidth',1);
% hold on;plot(sim_data_test.O,'--','LineWidth',1);
% f1=gcf;a1=gca;
% set(a1,'FontSize',14);
% xlabel('timestep [ ]');
% ylabel('[ ]');
% title(['Performance comparison - Net 1 (MSE: ', num2str(mse), ')']);
% legend('target output','system output')

figure;plot(Y_test,'r','LineWidth',1);
hold on;plot(sim_data_test.O,'--','LineWidth',1);
title(['Net 3 (MSE: ', num2str(Error), ')']);
xlabel('timestep');
legend('target output','system output');

disp(['MSE: ',num2str(Error)]);
%disp(['MSE HYS: ',num2str(mse_hys)]);
%d = sim_data;
%params = net.W;

save(['files/save_sim_data/Error_hill_N_' num2str(N)],'Error');
end