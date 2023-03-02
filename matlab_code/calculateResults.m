function [] = calculateResults(N)
% N is the amount of memory of the task defined
% e.g. N=2
% must run runTenNetworks first before calculating results


load('files/smallNetworks3.mat');
load(['files/smallNets_NARMA_hys/networks2_N_' num2str(N)]);
load(['files/smallNets_NARMA_hys/sim_data_N_' num2str(N)]);
load(['NARMA_data/NARMA' num2str(N)]);


% Data = learningTask(N,250000); % I/O supervised learning task

wash_out = 60000;
start = 80000-wash_out; 
len = 200000+start;
len_test = 15000;

% U = Data.u(start:len,1);  
% Y = Data.y(start:len,1);
U = dat.un(start:len,1);  
Y = dat.yn(start:len,1);

% prepare testing data
% U_test = Data.u(len+1:len+len_test,1);  
% Y_test = Data.y(len+1:len+len_test,1);
U_test = dat.un(len+1:len+len_test,1);  
Y_test = dat.yn(len+1:len+len_test,1);

Error = zeros(10,1);
net_fields = fieldnames(smallNetworks);
sim_data_fields = fieldnames(sim_data);
net2_fields = fieldnames(networks2);
NETS = [3 6];
for i =1:length(NETS)
    disp([' j = ' ,num2str(i) , ' of  10']);
    if (strcmp(smallNetworks.(net_fields{NETS(i)}).readout_type,'LENGTHS'))
        % sim_data.D lengths of every spring
        X = sim_data.(sim_data_fields{NETS(i)}).D(wash_out:end,:);  % throw first 100 steps away
    else
        % sim_data.Sx x-positions of each mass
        X = sim_data.(sim_data_fields{NETS(i)}).Sx(wash_out:end,:);  % throw first 100 steps away
    end

    Yw = Y(wash_out:end,:);

    W_out=X\Yw; % X*W_out=Yw e.g. (1000x77)*(77x1)=(1000x1) 


    net_test = networks2.(net2_fields{NETS(i)}); %output network
    net_test.W_out = W_out; % output weights for each spring


%     [net_test_out,sim_data_test] = simulate_ms_sys(net_test,U_test);
    [net_test_out,sim_data_test] = simulate_hysteresis_sys(net_test,U_test);
%     [net_test_out,sim_data_test] = simulate_maxwell_sys(net_test,U_test);
    Error(NETS(i)) = mean_squared_error(Y_test,sim_data_test.O);
    % disp(['MSE: ',num2str(error)])
end

save(['files/smallNets_NARMA_hys/Error_N_' num2str(N)],'Error');
%x=N*ones(10,1);
%figure; plot(x,Error,'x');

