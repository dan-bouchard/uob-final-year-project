function [networks] = make_ten_networks()

% load data
close all;
clear all;

data = init_ms_sys_data;
data.num = 24;
data.k_lim = [1 100;1 200];  
data.d_lim = [1 100;1 200];

networks.net1 = init_ms_sys_net(data);
networks.net2 = init_ms_sys_net(data);
networks.net3 = init_ms_sys_net(data);
networks.net4 = init_ms_sys_net(data);
networks.net5 = init_ms_sys_net(data);
networks.net6 = init_ms_sys_net(data);
networks.net7 = init_ms_sys_net(data);
networks.net8 = init_ms_sys_net(data);
networks.net9 = init_ms_sys_net(data);
networks.net10 = init_ms_sys_net(data);

