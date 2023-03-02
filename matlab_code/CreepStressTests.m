
k1 = 20;
k3 = 10;
d1 = 15;
d3 = 5;
kp1 = 7;
kp3 = 3;

time_step = 0.01;
end_time = 100;
end_time = 20;
t = linspace(0,end_time,end_time/time_step + 1);

linear_creep.x = zeros(size(t));
linear_creep.x_dot = zeros(size(t));
linear_creep.x2 = zeros(size(t));
linear_creep.x2_dot = zeros(size(t));
linear_creep.F = zeros(size(t));

linear_stress.x = zeros(size(t));
linear_stress.x_dot = zeros(size(t));
linear_stress.x2 = zeros(size(t));
linear_stress.x2_dot = zeros(size(t));
linear_stress.F = zeros(size(t));

nonlinear_creep.x = zeros(size(t));
nonlinear_creep.x_dot = zeros(size(t));
nonlinear_creep.x2 = zeros(size(t));
nonlinear_creep.x2_dot = zeros(size(t));
nonlinear_creep.F = zeros(size(t));

nonlinear_stress.x = zeros(size(t));
nonlinear_stress.x_dot = zeros(size(t));
nonlinear_stress.x2 = zeros(size(t));
nonlinear_stress.x2_dot = zeros(size(t));
nonlinear_stress.F = zeros(size(t));

dynamic_linear_creep.x = zeros(size(t));
dynamic_linear_creep.x_dot = zeros(size(t));
dynamic_linear_creep.x2 = zeros(size(t));
dynamic_linear_creep.x2_dot = zeros(size(t));
dynamic_linear_creep.F = zeros(size(t));

dynamic_linear_stress.x = zeros(size(t));
dynamic_linear_stress.x_dot = zeros(size(t));
dynamic_linear_stress.x2 = zeros(size(t));
dynamic_linear_stress.x2_dot = zeros(size(t));
dynamic_linear_stress.F = zeros(size(t));

dynamic_nonlinear_creep.x = zeros(size(t));
dynamic_nonlinear_creep.x_dot = zeros(size(t));
dynamic_nonlinear_creep.x2 = zeros(size(t));
dynamic_nonlinear_creep.x2_dot = zeros(size(t));
dynamic_nonlinear_creep.F = zeros(size(t));

dynamic_nonlinear_stress.x = zeros(size(t));
dynamic_nonlinear_stress.x_dot = zeros(size(t));
dynamic_nonlinear_stress.x2 = zeros(size(t));
dynamic_nonlinear_stress.x2_dot = zeros(size(t));
dynamic_nonlinear_stress.F = zeros(size(t));



p=5; B=1*0.05;
x_t = 0.03*ones(size(t));
Non_x_t = 0.03*ones(size(t));
strain_x_t = 3*ones(size(t));

f_t = B*sin(p*t);
strain_f_t = 500*B*sin(p*t);
df_t = B*p*cos(p*t);
strain_df_t = 500*B*p*cos(p*t);


dynamic_linear_stress.F(1) = k1*f_t(1);
dynamic_linear_stress.x_dot(1) = df_t(1);
dynamic_linear_stress.x(1) = f_t(1);

linear_stress.F(1) = k1*x_t(1);
%linear_stress.F(1) = k1*x_t(1) + kp1*x_t(1);
linear_stress.x_dot(1) = 0;
linear_stress.x(1) = x_t(1);
linear_stress.x2(1) = x_t(1);

nonlinear_stress.F(1) = k1*Non_x_t(1)+k3*Non_x_t(1)^3;
% nonlinear_stress.F(1) = k1*Non_x_t(1) + k3*Non_x_t(1)^3 + kp1*Non_x_t(1) + kp3*Non_x_t(1)^3;
nonlinear_stress.x_dot(1) = 0;
nonlinear_stress.x(1) = Non_x_t(1);
nonlinear_stress.x2(1) = Non_x_t(1);

dynamic_nonlinear_stress.F(1) = k1*f_t(1)+k3*f_t(1)^3;
dynamic_nonlinear_stress.x_dot(1) = df_t(1);
dynamic_nonlinear_stress.x(1) = f_t(1);
dynamic_nonlinear_stress.x2(1) = f_t(1);
for i=1:length(t)-1
    if(mod(i,1000)==0) % when to show steps
 	  disp([' i = ' ,num2str(i) , ' of ' , num2str(length(t)) ]);
    end
    %%% MAXWELL %%%
%     dynamic_linear_stress.F(i+1) = dynamic_linear_stress.F(i) + time_step*(k1*dynamic_linear_stress.x_dot(i)-(k1*dynamic_linear_stress.F(i))/d1);
%     dynamic_linear_stress.x_dot(i+1) = df_t(i+1);
%     dynamic_linear_stress.x(i+1) = f_t(i+1);
%     
%     linear_stress.F(i+1) = linear_stress.F(i) + time_step*(k1*linear_stress.x_dot(i)-(k1*linear_stress.F(i))/d1);
%     linear_stress.x(i+1) = x_t(i+1);
%     
%     linear_creep.F(i+1) = linear_creep.F(i) + time_step*(k1*linear_creep.x_dot(i)-(k1*linear_creep.F(i))/d1);
%     linear_creep.x_dot(i+1) = linear_creep.x_dot(i) + time_step*(strain_x_t(i)-linear_creep.F(i));
%     linear_creep.x(i+1) = linear_creep.x(i) + time_step*linear_creep.x_dot(i);
%     
%     dynamic_linear_creep.F(i+1) = dynamic_linear_creep.F(i) + time_step*(k1*dynamic_linear_creep.x_dot(i)-(k1*dynamic_linear_creep.F(i))/d1);
%     dynamic_linear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(strain_f_t(i)-dynamic_linear_creep.F(i));
%     dynamic_linear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
%     
%     nonlinear_creep.x_dot(i+1) = nonlinear_creep.x_dot(i) + time_step*(strain_x_t(i)-nonlinear_creep.F(i));
%     nonlinear_creep.x(i+1) = nonlinear_creep.x(i) + time_step*nonlinear_creep.x_dot(i);
%     max_fun = @(X2_dot)(d3*(nonlinear_creep.x_dot(i) - X2_dot)^3 + d1*(nonlinear_creep.x_dot(i) - X2_dot) - k3*nonlinear_creep.x2(i)^3 - k1*nonlinear_creep.x2(i));
%     nonlinear_creep.x2_dot(i+1) = fzero(max_fun, nonlinear_creep.x2_dot(i));
%     nonlinear_creep.x2(i+1) = nonlinear_creep.x2(i) + time_step*nonlinear_creep.x2_dot(i+1);
%     nonlinear_creep.F(i+1) = k3*nonlinear_creep.x2(i+1)^3+k1*nonlinear_creep.x2(i+1);
%     
%      
%     dynamic_nonlinear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(f_t(i)-dynamic_linear_creep.F(i));
%     dynamic_nonlinear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
%     max_fun = @(X2_dot)(d3*(dynamic_nonlinear_creep.x_dot(i) - X2_dot)^3 + d1*(dynamic_nonlinear_creep.x_dot(i) - X2_dot) - k3*dynamic_nonlinear_creep.x2(i)^3 - k1*dynamic_nonlinear_creep.x2(i));
%     dynamic_nonlinear_creep.x2_dot(i+1) = fzero(max_fun, dynamic_nonlinear_creep.x2_dot(i));
%     dynamic_nonlinear_creep.x2(i+1) = dynamic_nonlinear_creep.x2(i) + time_step*dynamic_nonlinear_creep.x2_dot(i+1);
%     dynamic_nonlinear_creep.F(i+1) = k3*dynamic_nonlinear_creep.x2(i+1)^3+k1*dynamic_nonlinear_creep.x2(i+1);
%     
%     nonlinear_stress.x(i+1) = Non_x_t(i+1);
%     max_fun = @(X2_dot)(d3*(nonlinear_stress.x_dot(i) - X2_dot)^3 + d1*(nonlinear_stress.x_dot(i) - X2_dot) - k3*nonlinear_stress.x2(i)^3 - k1*nonlinear_stress.x2(i));
%     nonlinear_stress.x2_dot(i+1) = fzero(max_fun, nonlinear_stress.x2_dot(i));
%     nonlinear_stress.x2(i+1) = nonlinear_stress.x2(i) + time_step*nonlinear_stress.x2_dot(i+1);
%     nonlinear_stress.F(i+1) = k3*nonlinear_stress.x2(i+1)^3+k1*nonlinear_stress.x2(i+1);
%     
%     dynamic_nonlinear_stress.x(i+1) = f_t(i+1);
%     dynamic_nonlinear_stress.x_dot(i+1) = df_t(i+1);
%     max_fun = @(X2_dot)(d3*(dynamic_nonlinear_stress.x_dot(i) - X2_dot)^3 + d1*(dynamic_nonlinear_stress.x_dot(i) - X2_dot) - k3*dynamic_nonlinear_stress.x2(i)^3 - k1*dynamic_nonlinear_stress.x2(i));
%     dynamic_nonlinear_stress.x2_dot(i+1) = fzero(max_fun, dynamic_nonlinear_stress.x2_dot(i));
%     dynamic_nonlinear_stress.x2(i+1) = dynamic_nonlinear_stress.x2(i) + time_step*dynamic_nonlinear_stress.x2_dot(i+1);
%     dynamic_nonlinear_stress.F(i+1) = k3*dynamic_nonlinear_stress.x2(i+1)^3+k1*dynamic_nonlinear_stress.x2(i+1);
    
    %%%% VOIGT %%%%    
    dynamic_linear_stress.x_dot(i+1) = df_t(i+1);
    dynamic_linear_stress.x(i+1) = f_t(i+1);
    dynamic_linear_stress.F(i+1) = k1*dynamic_linear_stress.x(i+1) + d1*dynamic_linear_stress.x_dot(i+1);
    
    linear_stress.x(i+1) = x_t(i+1);
    linear_stress.F(i+1) = k1*linear_stress.x(i+1);
        
    linear_creep.x_dot(i+1) = linear_creep.x_dot(i) + time_step*(strain_x_t(i)-linear_creep.F(i));
    linear_creep.x(i+1) = linear_creep.x(i) + time_step*linear_creep.x_dot(i);
    linear_creep.F(i+1) = k1*linear_creep.x(i+1) + d1*linear_creep.x_dot(i+1);    
    
    dynamic_linear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(strain_f_t(i)-dynamic_linear_creep.F(i));
    dynamic_linear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
    dynamic_linear_creep.F(i+1) = k1*dynamic_linear_creep.x(i+1) + d1*dynamic_linear_creep.x_dot(i+1);
    
    nonlinear_creep.x_dot(i+1) = nonlinear_creep.x_dot(i) + time_step*(strain_x_t(i)-nonlinear_creep.F(i));
    nonlinear_creep.x(i+1) = nonlinear_creep.x(i) + time_step*nonlinear_creep.x_dot(i);
    nonlinear_creep.F(i+1) = k1*nonlinear_creep.x(i+1) + k3*nonlinear_creep.x(i+1)^3 + d1*nonlinear_creep.x_dot(i+1) + d3*nonlinear_creep.x_dot(i+1)^3;
     
    dynamic_nonlinear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(f_t(i)-dynamic_linear_creep.F(i));
    dynamic_nonlinear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
    dynamic_nonlinear_creep.F(i+1) = k1*dynamic_nonlinear_creep.x(i+1) + k3*dynamic_nonlinear_creep.x(i+1)^3 + d1*dynamic_nonlinear_creep.x_dot(i+1) + d3*dynamic_nonlinear_creep.x_dot(i+1)^3;
    
    nonlinear_stress.x(i+1) = Non_x_t(i+1);
    nonlinear_stress.F(i+1) = k1*nonlinear_stress.x(i+1) + k3*nonlinear_stress.x(i+1).^3 + d1*nonlinear_stress.x_dot(i+1) + d3*nonlinear_stress.x_dot(i+1).^3;
    
    dynamic_nonlinear_stress.x(i+1) = f_t(i+1);
    dynamic_nonlinear_stress.x_dot(i+1) = df_t(i+1);
    dynamic_nonlinear_stress.F(i+1) = k1*dynamic_nonlinear_stress.x(i+1) + k3*dynamic_nonlinear_stress.x(i+1)^3 + d1*dynamic_nonlinear_stress.x_dot(i+1) + d3*dynamic_nonlinear_stress.x_dot(i+1)^3;
        
    %%% KELVIN_A %%%    
    
%     dynamic_linear_stress.x_dot(i+1) = df_t(i+1);
%     dynamic_linear_stress.x(i+1) = f_t(i+1);
%     dynamic_linear_stress.x2_dot(i+1) = dynamic_linear_stress.x_dot(i) - (k1/d1)*dynamic_linear_stress.x2(i);
%     dynamic_linear_stress.x2(i+1) = dynamic_linear_stress.x2(i) + time_step*dynamic_linear_stress.x2_dot(i+1);
%     dynamic_linear_stress.F(i+1) = k1*dynamic_linear_stress.x2(i+1) + kp1*dynamic_linear_stress.x(i+1);
%     
%     linear_stress.x(i+1) = x_t(i+1);
%     linear_stress.x2_dot(i+1) = linear_stress.x_dot(i) - (k1/d1)*linear_stress.x2(i);
%     linear_stress.x2(i+1) = linear_stress.x2(i) + time_step*linear_stress.x2_dot(i+1);
%     linear_stress.F(i+1) = k1*linear_stress.x2(i+1) + kp1*linear_stress.x(i+1);
%     
%     linear_creep.x_dot(i+1) = linear_creep.x_dot(i) + time_step*(strain_x_t(i)-linear_creep.F(i));
%     linear_creep.x(i+1) = linear_creep.x(i) + time_step*linear_creep.x_dot(i);
%     linear_creep.x2_dot(i+1) = linear_creep.x_dot(i) - (k1/d1)*linear_creep.x2(i);
%     linear_creep.x2(i+1) = linear_creep.x2(i) + time_step*linear_creep.x2_dot(i+1);
%     linear_creep.F(i+1) = k1*linear_creep.x2(i+1) + kp1*linear_creep.x(i+1);
%     
%     dynamic_linear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(strain_f_t(i)-dynamic_linear_creep.F(i));
%     dynamic_linear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
%     dynamic_linear_creep.x2_dot(i+1) = dynamic_linear_creep.x_dot(i) - (k1/d1)*dynamic_linear_creep.x2(i);
%     dynamic_linear_creep.x2(i+1) = dynamic_linear_creep.x2(i) + time_step*dynamic_linear_creep.x2_dot(i+1);
%     dynamic_linear_creep.F(i+1) = k1*dynamic_linear_creep.x2(i+1) + kp1*dynamic_linear_creep.x(i+1);
%     
%     nonlinear_creep.x_dot(i+1) = nonlinear_creep.x_dot(i) + time_step*(strain_x_t(i)-nonlinear_creep.F(i));
%     nonlinear_creep.x(i+1) = nonlinear_creep.x(i) + time_step*nonlinear_creep.x_dot(i);
%     kelvina_fun = @(X2_dot)(d3*(nonlinear_creep.x_dot(i) - X2_dot)^3 + d1*(nonlinear_creep.x_dot(i) - X2_dot) - k3*nonlinear_creep.x2(i)^3 - k1*nonlinear_creep.x2(i));
%     nonlinear_creep.x2_dot(i+1) = fzero(kelvina_fun, nonlinear_creep.x2_dot(i));
%     nonlinear_creep.x2(i+1) = nonlinear_creep.x2(i) + time_step*nonlinear_creep.x2_dot(i+1);
%     nonlinear_creep.F(i+1) = k1*nonlinear_creep.x2(i+1) + k3*nonlinear_creep.x2(i+1)^3 + kp1*nonlinear_creep.x(i+1) + kp3*nonlinear_creep.x(i+1)^3;
%      
%     dynamic_nonlinear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(f_t(i)-dynamic_linear_creep.F(i));
%     dynamic_nonlinear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
%     kelvina_fun = @(X2_dot)(d3*(dynamic_nonlinear_creep.x_dot(i) - X2_dot)^3 + d1*(dynamic_nonlinear_creep.x_dot(i) - X2_dot) - k3*dynamic_nonlinear_creep.x2(i)^3 - k1*dynamic_nonlinear_creep.x2(i));
%     dynamic_nonlinear_creep.x2_dot(i+1) = fzero(kelvina_fun, dynamic_nonlinear_creep.x2_dot(i));
%     dynamic_nonlinear_creep.x2(i+1) = dynamic_nonlinear_creep.x2(i) + time_step*dynamic_nonlinear_creep.x2_dot(i+1);
%     dynamic_nonlinear_creep.F(i+1) = k1*dynamic_nonlinear_creep.x2(i+1) + k3*dynamic_nonlinear_creep.x2(i+1)^3 + kp1*dynamic_nonlinear_creep.x(i+1) + kp3*dynamic_nonlinear_creep.x(i+1)^3;
%     
%     nonlinear_stress.x(i+1) = Non_x_t(i+1);
%     kelvina_fun = @(X2_dot)(d3*(nonlinear_stress.x_dot(i) - X2_dot)^3 + d1*(nonlinear_stress.x_dot(i) - X2_dot) - k3*nonlinear_stress.x2(i)^3 - k1*nonlinear_stress.x2(i));
%     nonlinear_stress.x2_dot(i+1) = fzero(kelvina_fun, nonlinear_stress.x2_dot(i));
%     nonlinear_stress.x2(i+1) = nonlinear_stress.x2(i) + time_step*nonlinear_stress.x2_dot(i+1);
%     nonlinear_stress.F(i+1) = k1*nonlinear_stress.x2(i+1) + k3*nonlinear_stress.x2(i+1)^3 + kp1*nonlinear_stress.x(i+1) + kp3*nonlinear_stress.x(i+1)^3; 
%     
%     dynamic_nonlinear_stress.x(i+1) = f_t(i+1);
%     dynamic_nonlinear_stress.x_dot(i+1) = df_t(i+1);
%     kelvina_fun = @(X2_dot)(d3*(dynamic_nonlinear_stress.x_dot(i) - X2_dot)^3 + d1*(dynamic_nonlinear_stress.x_dot(i) - X2_dot) - k3*dynamic_nonlinear_stress.x2(i)^3 - k1*dynamic_nonlinear_stress.x2(i));
%     dynamic_nonlinear_stress.x2_dot(i+1) = fzero(kelvina_fun, dynamic_nonlinear_stress.x2_dot(i));
%     dynamic_nonlinear_stress.x2(i+1) = dynamic_nonlinear_stress.x2(i) + time_step*dynamic_nonlinear_stress.x2_dot(i+1);
%     dynamic_nonlinear_stress.F(i+1) = k1*dynamic_nonlinear_stress.x2(i+1) + k3*dynamic_nonlinear_stress.x2(i+1)^3 + kp1*dynamic_nonlinear_stress.x(i+1) + kp3*dynamic_nonlinear_stress.x(i+1)^3;
    
    %%%%% KELVIN B %%%%%

%     dynamic_linear_stress.x_dot(i+1) = df_t(i+1);
%     dynamic_linear_stress.x(i+1) = f_t(i+1);
%     dynamic_linear_stress.x2_dot(i+1) = dynamic_linear_stress.x_dot(i) - (k1/d1)*dynamic_linear_stress.x2(i) + (kp1/d1)*(dynamic_linear_stress.x(i) - dynamic_linear_stress.x2(i));
%     dynamic_linear_stress.x2(i+1) = dynamic_linear_stress.x2(i) + time_step*dynamic_linear_stress.x2_dot(i+1);
%     dynamic_linear_stress.F(i+1) = k1*dynamic_linear_stress.x2(i+1);
%     
%     linear_stress.x(i+1) = x_t(i+1);
%     linear_stress.x2_dot(i+1) = linear_stress.x_dot(i) - (k1/d1)*linear_stress.x2(i) + (kp1/d1)*(linear_stress.x(i) - linear_stress.x2(i));
%     linear_stress.x2(i+1) = linear_stress.x2(i) + time_step*linear_stress.x2_dot(i+1);
%     linear_stress.F(i+1) = k1*linear_stress.x2(i+1);
%     
%     linear_creep.x_dot(i+1) = linear_creep.x_dot(i) + time_step*(strain_x_t(i)-linear_creep.F(i));
%     linear_creep.x(i+1) = linear_creep.x(i) + time_step*linear_creep.x_dot(i);
%     linear_creep.x2_dot(i+1) = linear_creep.x_dot(i) - (k1/d1)*linear_creep.x2(i) + (kp1/d1)*(linear_creep.x(i) - linear_creep.x2(i));
%     linear_creep.x2(i+1) = linear_creep.x2(i) + time_step*linear_creep.x2_dot(i+1);
%     linear_creep.F(i+1) = k1*linear_creep.x2(i+1);
%     
%     dynamic_linear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(strain_f_t(i)-dynamic_linear_creep.F(i));
%     dynamic_linear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
%     dynamic_linear_creep.x2_dot(i+1) = dynamic_linear_creep.x_dot(i) - (k1/d1)*dynamic_linear_creep.x2(i) + (kp1/d1)*(dynamic_linear_creep.x(i) - dynamic_linear_creep.x2(i));
%     dynamic_linear_creep.x2(i+1) = dynamic_linear_creep.x2(i) + time_step*dynamic_linear_creep.x2_dot(i+1);
%     dynamic_linear_creep.F(i+1) = k1*dynamic_linear_creep.x2(i+1);
%     
%     nonlinear_creep.x_dot(i+1) = nonlinear_creep.x_dot(i) + time_step*(strain_x_t(i)-nonlinear_creep.F(i));
%     nonlinear_creep.x(i+1) = nonlinear_creep.x(i) + time_step*nonlinear_creep.x_dot(i);
%     kelvinb_fun = @(X2_dot)(d3*(nonlinear_creep.x_dot(i) - X2_dot)^3 + d1*(nonlinear_creep.x_dot(i) - X2_dot) + kp3*(nonlinear_creep.x(i) - nonlinear_creep.x2(i))^3 + kp1*(nonlinear_creep.x(i) - nonlinear_creep.x2(i)) - k3*nonlinear_creep.x2(i)^3 - k1*nonlinear_creep.x2(i));
%     nonlinear_creep.x2_dot(i+1) = fzero(kelvinb_fun, nonlinear_creep.x2_dot(i));
%     nonlinear_creep.x2(i+1) = nonlinear_creep.x2(i) + time_step*nonlinear_creep.x2_dot(i+1);
%     nonlinear_creep.F(i+1) = k1*nonlinear_creep.x2(i+1) + k3*nonlinear_creep.x2(i+1)^3;
%      
%     dynamic_nonlinear_creep.x_dot(i+1) = dynamic_linear_creep.x_dot(i) + time_step*(f_t(i)-dynamic_linear_creep.F(i));
%     dynamic_nonlinear_creep.x(i+1) = dynamic_linear_creep.x(i) + time_step*dynamic_linear_creep.x_dot(i);
%     kelvinb_fun = @(X2_dot)(d3*(dynamic_linear_creep.x_dot(i) - X2_dot)^3 + d1*(dynamic_linear_creep.x_dot(i) - X2_dot) + kp3*(dynamic_linear_creep.x(i) - dynamic_linear_creep.x2(i))^3 + kp1*(dynamic_linear_creep.x(i) - dynamic_linear_creep.x2(i)) - k3*dynamic_linear_creep.x2(i)^3 - k1*dynamic_linear_creep.x2(i));
%     dynamic_nonlinear_creep.x2_dot(i+1) = fzero(kelvinb_fun, dynamic_nonlinear_creep.x2_dot(i));
%     dynamic_nonlinear_creep.x2(i+1) = dynamic_nonlinear_creep.x2(i) + time_step*dynamic_nonlinear_creep.x2_dot(i+1);
%     dynamic_nonlinear_creep.F(i+1) = k1*dynamic_nonlinear_creep.x2(i+1) + k3*dynamic_nonlinear_creep.x2(i+1)^3;
%     
%     nonlinear_stress.x(i+1) = Non_x_t(i+1);
%     kelvinb_fun = @(X2_dot)(d3*(nonlinear_stress.x_dot(i) - X2_dot)^3 + d1*(nonlinear_stress.x_dot(i) - X2_dot) + kp3*(nonlinear_stress.x(i) - nonlinear_stress.x2(i))^3 + kp1*(nonlinear_stress.x(i) - nonlinear_stress.x2(i)) - k3*nonlinear_stress.x2(i)^3 - k1*nonlinear_stress.x2(i));
%     nonlinear_stress.x2_dot(i+1) = fzero(kelvinb_fun, nonlinear_stress.x2_dot(i));
%     nonlinear_stress.x2(i+1) = nonlinear_stress.x2(i) + time_step*nonlinear_stress.x2_dot(i+1);
%     nonlinear_stress.F(i+1) = k1*nonlinear_stress.x2(i+1) + k3*nonlinear_stress.x2(i+1)^3; 
%     
%     dynamic_nonlinear_stress.x(i+1) = f_t(i+1);
%     dynamic_nonlinear_stress.x_dot(i+1) = df_t(i+1);
%     kelvinb_fun = @(X2_dot)(d3*(dynamic_nonlinear_stress.x_dot(i) - X2_dot)^3 + d1*(dynamic_nonlinear_stress.x_dot(i) - X2_dot) + kp3*(dynamic_nonlinear_stress.x(i) - dynamic_nonlinear_stress.x2(i))^3 + kp1*(dynamic_nonlinear_stress.x(i) - dynamic_nonlinear_stress.x2(i)) - k3*dynamic_nonlinear_stress.x2(i)^3 - k1*dynamic_nonlinear_stress.x2(i));
%     dynamic_nonlinear_stress.x2_dot(i+1) = fzero(kelvinb_fun, dynamic_nonlinear_stress.x2_dot(i));
%     dynamic_nonlinear_stress.x2(i+1) = dynamic_nonlinear_stress.x2(i) + time_step*dynamic_nonlinear_stress.x2_dot(i+1);
%     dynamic_nonlinear_stress.F(i+1) = k1*dynamic_nonlinear_stress.x2(i+1) + k3*dynamic_nonlinear_stress.x2(i+1)^3;
    
end

figure; plot(t,linear_stress.F); hold on;
plot(t,linear_creep.x); 
%plot(t,nonlinear_stress.F);
plot(t,nonlinear_creep.x); 
hold off;
legend('Stress','Strain','Nonlinear Strain');
%'Nonlinear Stress','Nonlinear Strain');
xlabel('Time (s)'); ylabel('Stress (Pa), Strain'); title('Creep and Stress Relaxation Tests for Kelvin B Model');

figure; plot(t,dynamic_linear_stress.F); hold on;
plot(t,dynamic_linear_creep.x);
%plot(t,dynamic_nonlinear_stress.F);
%plot(t,dynamic_nonlinear_creep.x); 
hold off;
legend('Stress','Strain');
%,'Nonlinear Stress','Nonlinear Strain');
xlabel('Time (s)'); ylabel('Stress (Pa), Strain'); title('Dynamic stress and strain responses for Kelvin  Model');