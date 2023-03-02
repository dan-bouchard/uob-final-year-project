% learning simple memory task

% load data
close all;
clear all;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  making net 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = init_ms_sys_data;

% change parameters from default
% to ones that are appropiate for this task
data.num = 20; 		
data.show_steps = 1000;

% define ranges for the randomly intialized
% dynamic parameters of the springs stiffness and damping
data.k_lim = [1 100;1 200];  
data.d_lim = [1 100;1 200];


data.show_plot = 1;
%  data.readout_type = 'POSITIONS';
data.readout_type = 'LENGTHS';


% initialize a random net with given values
net = init_ms_sys_net(data); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loading data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data = learningTask(2,250000); % I/O supervised learning task

wash_out = 60000;
start = 80000-wash_out; 
len = 200000+start;
len_test = 15000;

U = Data.u(start:len,1);  
Y = Data.y(start:len,1);

% prepare testing data
U_test = Data.u(len+1:len+len_test,1);  
Y_test = Data.y(len+1:len+len_test,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  simulating net 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[net2,sim_data] = simulate_ms_sys(net,U);
 
% learn output weights with linear regression
if (strcmp(net.readout_type,'LENGTHS'))
    % sim_data.D lengths of every spring
	X = sim_data.D(wash_out:end,:);  % throw first 100 steps away
else
    % sim_data.Sx x-positions of each mass
	X = sim_data.Sx(wash_out:end,:);  % throw first 100 steps away
end

Yw = Y(wash_out:end,:);

W_out=X\Yw; % X*W_out=Yw e.g. (1000x77)*(77x1)=(1000x1) 


net_test = net2; %output network
net_test.W_out = W_out; % output weights for each spring
o = X*W_out;

% in case you want to test how good the weights 
% represent the learning data
% figure;plot(o); hold on;plot(Yw,'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% testing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[net_test_out,sim_data_test] = simulate_ms_sys(net_test,U_test);



% plot results
figure;plot(Y_test,'r','LineWidth',1);
hold on;plot((mapstd(sim_data_test.O'))','--','LineWidth',1);
f1=gcf;a1=gca;
set(a1,'FontSize',14);
xlabel('timestep [ ]');
ylabel('[ ]');
title('Performance comparison')
legend('target output','system output')


% caculate and print MSE
disp(['MSE: ',num2str(mean_squared_error(Y_test,(mapstd(sim_data_test.O'))'))])


% plot the structure of the used network
plot_graph(net_test)
