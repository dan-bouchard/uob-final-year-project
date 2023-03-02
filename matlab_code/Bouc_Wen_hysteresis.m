function [F, z_new] = Bouc_Wen_hysteresis(time_step, x, x_dot, z, beta, gamma, k3, k1, d3, d1)

% Vectors setup
% x = zeros(size(t));
% y = zeros(size(t)); % x dot
% z = zeros(size(t));
% v = zeros(size(t));

% parameters
% p=1; B=1;
A=1;
alpha = 0; n = 1.2;
% doesn't seem to work with beta<0
% beta = 0.25; gamma = -0.75;

   
% x_2 = rk4(@(x)( y ), x, time_step);
% y_2 = rk4(@(y)( -d*y-alpha*k*x-(1-alpha)*k*z ), y, time_step);
% z_2 = rk4(@(z)( A*y-beta*abs(y)*abs(z).^(n-1)*z-gamma*y*abs(z).^n ) , z, time_step);    


%transient = 30;
F_s = alpha*(k3*x^3+k1*x) + d3*x_dot^3+d1*x_dot;
F_z = (1-alpha)*(k3*z^3+k1*z);
F = F_s + F_z;

z_new = z + time_step*( A*x_dot-beta*abs(x_dot)*abs(z).^(n-1)*z-gamma*x_dot*abs(z).^n );
% z_new = rk4(@(z)( A*x_dot-beta*abs(x_dot)*abs(z).^(n-1)*z-gamma*x_dot*abs(z).^n ) , z, time_step);   

% F = alpha*k*x_2+(1-alpha)*k*z_2+d*y_2;
%plot(x(transient*1/time_step:end),z(transient*1/time_step:end)); 
%xlabel('x'); ylabel('z'); title(['Bouc-Wen Hysteresis Loop (\beta=' num2str(beta) ', \gamma=' num2str(gamma) ')']);
%grid 'on'


