N = [5 10 15 25 40 60 100];
% n = [5 10 15 20];
% N = [2 3 4 5];
% N = [5 10 15 20];
N_minus = N-2; %N_minus = N-0.5; N_minus = N-0.1;
N_plus = N+2; %N_plus = N+0.5; N_plus = N+0.1;
N_Plus = N+3.5;

colour_1 = [27/255 200/255 119/255]; % murky green
colour_2 = [240/255 95/255 2/255]; % orange/red
colour_3 = [110/255 100/255 250/255]; % blue/purple
colour_4 = [231/255 41/255 138/255]; % dark pink
colour_5 = [102/255 166/255 30/255]; % light green
colour_6 = [166/255 118/255 29/255]; % brown
colour_7 = [102/255 102/255 102/255]; % grey

Col_1 = [228/255 26/255 28/255];
Col_2 = [55/255 126/255 184/255];
Col_3 = [77/255 175/255 74/255];
Col_4 = [152/255 78/255 163/255];
Col_5 = [255/255 127/255 0];
Col_6 = [166/255 86/255 40/255];
Col_7 = [231/255 41/255 138/255];

% load('Data\hill_activation_memory'); hill_activation = hill_activation_memory;
load('Data\hill_memory'); hill = hill_memory;
load('Data\hys_memory'); hys = hys_memory;
load('Data\kelvin_a_memory'); kelvin_a = kelvin_a_memory;
load('Data\kelvin_b_memory'); kelvin_b = kelvin_b_memory;
load('Data\max_memory'); max = max_memory;
load('Data\voigt_memory'); voigt = voigt_memory;
load('Data/hill_small_CE_memory'); hill_small_CE = hill_small_CE_memory;

load('Data/kelvin_a_linear_memory'); kelvin_a_linear = kelvin_a_linear_memory;
load('Data/kelvin_b_linear_memory'); kelvin_b_linear = kelvin_b_linear_memory;
load('Data/max_linear_memory'); max_linear = max_linear_memory;

load('Data/hill_small_CE_memory_N_5');
load('Data/hill_Two_Thirds_memory_N_5');
load('Data/hill_Full_Input_memory_N_5');


x = ones(10,length(N));
for i=1:length(N)
    x(:,i) = x(:,i)*N(i);
end
x = reshape(x,10*length(N),1);

hill_mean = mean(hill);
%hill_activation_mean = mean(hill_activation);
hys_mean = mean(hys);
kelvin_a_mean = mean(kelvin_a);
kelvin_b_mean = mean(kelvin_b);
max_mean = mean(max);
voigt_mean = mean(voigt);
kelvin_a_linear_mean = mean(kelvin_a_linear);
kelvin_b_linear_mean = mean(kelvin_b_linear);
max_linear_mean = mean(max_linear);
hill_small_CE_mean = mean(hill_small_CE);
%hill_small_CE_mean = mean(hill_small_CE_memory_N_5);
hill_Two_Thirds_mean = mean(hill_Two_Thirds_memory_N_5);
hill_Full_Input_mean = mean(hill_Full_Input_memory_N_5);

hill_std = std(hill);
%hill_activation_std = std(hill_activation);
hys_std = std(hys);
kelvin_a_std = std(kelvin_a);
kelvin_b_std = std(kelvin_b);
max_std = std(max);
voigt_std = std(voigt);
kelvin_a_linear_std = std(kelvin_a_linear);
kelvin_b_linear_std = std(kelvin_b_linear);
max_linear_std = std(max_linear);
hill_small_CE_std = std(hill_small_CE);

%hill_small_CE_std = std(hill_small_CE_memory_N_5);
hill_Two_Thirds_std = std(hill_Two_Thirds_memory_N_5);
hill_Full_Input_std = std(hill_Full_Input_memory_N_5);

hill = reshape(hill,10*length(N),1);
% hill_activation = reshape(hill_activation,10*length(N),1);
hys = reshape(hys,10*length(N),1);
kelvin_a = reshape(kelvin_a,10*length(N),1);
kelvin_b = reshape(kelvin_b,10*length(N),1);
max = reshape(max,10*length(N),1);
voigt = reshape(voigt,10*length(N),1);
kelvin_a_linear = reshape(kelvin_a_linear,10*length(N),1);
kelvin_b_linear = reshape(kelvin_b_linear,10*length(N),1);
max_linear = reshape(max_linear,10*length(N),1);
hill_small_CE = reshape(hill_small_CE,10*length(N),1);

%Data = [voigt max kelvin_a kelvin_b hill hill_activation hys kelvin_a_linear kelvin_b_linear max_linear];
%Data = [voigt max kelvin_a kelvin_b hill hys];
Data = [voigt max kelvin_a kelvin_b hill hys kelvin_a_linear kelvin_b_linear max_linear hill_small_CE];
%Mean = [voigt_mean; max_mean; kelvin_a_mean; kelvin_b_mean; hill_mean; hill_activation_mean; hys_mean; kelvin_a_linear_mean; kelvin_b_linear_mean; max_linear_mean];
%Mean = [voigt_mean; max_mean; kelvin_a_mean; kelvin_b_mean; hill_mean; hys_mean];
Mean = [voigt_mean; max_mean; kelvin_a_mean; kelvin_b_mean; hill_mean; hys_mean; kelvin_a_linear_mean; kelvin_b_linear_mean; max_linear_mean; hill_small_CE_mean];
%Std =  [voigt_std; max_std; kelvin_a_std; kelvin_b_std; hill_std; hill_activation_std; hys_std; kelvin_a_linear_std; kelvin_b_linear_std; max_linear_std];
%Std =  [voigt_std; max_std; kelvin_a_std; kelvin_b_std; hill_std; hys_std];
Std =  [voigt_std; max_std; kelvin_a_std; kelvin_b_std; hill_std; hys_std; kelvin_a_linear_std; kelvin_b_linear_std; max_linear_std; hill_small_CE_std];

idx = [5 10];
figure; %plot(x,Data(:,idx),'x', 'MarkerSize',8);
plot(x,Data(:,idx(1)),'rx', 'MarkerSize',8); hold on; plot(x,Data(:,idx(2)),'b*', 'MarkerSize',8);
%plot(n(1)*ones(10,1),Data(1:10,5),'rx', 'MarkerSize',8); hold on;
%plot(n(2)*ones(10,1),hill_small_CE_memory_N_5,'bx', 'MarkerSize',8);
%plot(n(3)*ones(10,1),hill_Two_Thirds_memory_N_5,'gx', 'MarkerSize',8);
%plot(n(4)*ones(10,1),hill_Full_Input_memory_N_5,'mx', 'MarkerSize',8);
%plot(x,Data(:,idx(3)),'go', 'MarkerSize',8);

errorbar(N_minus, Mean(idx(1),:), Std(idx(1),:),'LineStyle','none', 'Marker','.','MarkerSize',16,'Color','red'); hold on;
errorbar(N_plus, Mean(idx(2),:), Std(idx(2),:),'LineStyle','none','Marker','.','MarkerSize',16,'Color','blue');
%errorbar(N_Plus, Mean(idx(3),:), Std(idx(3),:),'LineStyle','none','Marker','.','MarkerSize',16,'Color','green');
hold off;

xlabel('N'); ylabel('MSE'); xlim([0 110]); %ylim([1e-6 1]); title('Memory Task Results');
%xlim([0 25]); ylim([0.01 2]); title('Memory Task Results');
%xlim([0 25]); title('NARMA Task Results');
title('Memory Task Results');

a1 = gca;
set(a1, 'YScale', 'log');
set(a1,'FontSize',12);
legend('Hill', 'Hill Small CE');
%legend('Hill','Hill with small CE','Two Thirds Input Connections','Full Input Connections');
%legend('Voigt', 'Maxwell', 'Kelvin A', 'Kelvin B', 'Hill', 'Hill with activation', 'Hysteresis');
