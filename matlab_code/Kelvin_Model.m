% Kelvin model

time_step = 0.01;
end_time = 100;
t = linspace(0,end_time,end_time/time_step + 1);
x = zeros(size(t));
x_dot = zeros(size(t));
x2 = zeros(size(t));
x2_dot = zeros(size(t));
x1_dot = zeros(size(t));
F_see = zeros(size(t));
F_pee = zeros(size(t));
F = zeros(size(t));

p=1; B=1; 
ks=20; d=10; kp=20;
% ks1 = 20; ks3 = 10;
% kp1 = 2; kp3 = 1;
% d1 = 10; d3 = 4;

% sinusoidal forcing
f_t = B*sin(p*t);
% step forcing
f_t = ones(size(t));

% Euler Integration

for i=1:length(t)-1
%     v2(i+1) = v2(i) + time_step*(f_t(i)-(ks3*(x2(i)-x1(i)).^3 + ks1*(x2(i)-x1(i)) + kp3*(x2(i)).^3 + kp1*(x2(i)) ));
%     x2(i+1) = x2(i) + v2(i)*time_step;
%     v1(i+1) = v1(i) + time_step*(ks3*(x2(i)-x1(i)).^3 + ks1*(x2(i)-x1(i)) - d*v1(i));
%     x1(i+1) = x1(i) + v1(i)*time_step;   
    
    x_dot(i+1) = x_dot(i) + time_step*(f_t(i)-F(i));
    x(i+1) = x(i) + time_step*x_dot(i);
    
    x2_dot(i+1) = x_dot(i) - ks*x2(i)/d;
%     fun = @(X2_dot)(d3*(x_dot(i) - X2_dot)^3 + d1*(x_dot(i) - X2_dot) - k3*x2(i)^3 - k1*x2(i));
%     x2_dot(i+1) = fzero(fun, x2_dot(i));
    
    x2(i+1) = x2(i) + time_step*x2_dot(i);
    x1_dot(i+1) = x_dot(i+1)-x2_dot(i+1);
    
    Fs = ks*x2(i+1);
%     Fs = ks1*x2(i+1)+ks3*x2(i+1)^3;
    
    Fd = d*(x_dot(i+1) - x2_dot(i+1));
%     Fd = d1*(x_dot(i+1) - x2_dot(i+1)) + d3*(x_dot(i+1) - x2_dot(i+1))^3;
    
    F_see(i+1) = 0.5*(Fs+Fd);
    
    F_pee(i+1) = kp*x(i+1);
%     F_pee(i+1) = kp1*x(i+1) + kp3*x(i+1)^3;
    F(i+1) = F_see(i+1)+F_pee(i+1);
end

% Plot Force v Extension

figure; plot(x(20/time_step:end),F(20/time_step:end)); xlabel('Extension (m)'); ylabel('Force (N)'); title('Nonlinear Kelvin Model A');
figure; plot(x,F); xlabel('Extension (m)'); ylabel('Force (N)'); title('Nonlinear Kelvin Model A');


