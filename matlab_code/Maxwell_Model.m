% Maxwell model

time_step = 0.01;
end_time = 100;
t = linspace(0,end_time,end_time/time_step + 1);
% x1 = zeros(size(t));
x2 = zeros(size(t));
% v1 = zeros(size(t));
% v2 = zeros(size(t));
x2_dot = zeros(size(t));
F = zeros(size(t));
x = zeros(size(t));
x_dot = zeros(size(t));


p=2; B=1; k=5; d=10;
% k1 = 20; k3 = 10;
% d1 = 10; d3 = 7;

k1 = 20; k3 = 10;
d1 = 10; d3 = 20;



% sinusoidal forcing
f_t = B*sin(p*t);
% step forcing
f_t = ones(size(t));

% Euler Integration

for i=1:length(t)-1
%     v2(i+1) = v2(i) + time_step*(f_t(i)-(k3*(x2(i)-x1(i)).^3 + k1*(x2(i)-x1(i)) ));
%     x2(i+1) = x2(i) + v2(i)*time_step;
%     v1(i+1) = v1(i) + time_step*(k3*(x2(i)-x1(i)).^3 + k1*(x2(i)-x1(i)) - d*v1(i));
%     x1(i+1) = x1(i) + v1(i)*time_step;

%     F(i+1) = F(i) + time_step*(k*x_dot(i)-k*F(i)/d);
%     x_dot(i+1) = x_dot(i) + time_step*(f_t(i)-F(i));
%     x(i+1) = x(i) + time_step*x_dot(i);
    
%     x2_dot(i)=x_dot(i)-k*x2(i)/d;
%     F(i)=k*x2(i);
%     x_dot(i+1) = x_dot(i) + time_step*(f_t(i)-F(i));
%     x(i+1) = x(i) + time_step*x_dot(i);
%     x2(i+1) = x2(i) + time_step*x2_dot(i);
    
    x_dot(i+1) = x_dot(i) + time_step*(f_t(i)-F(i));
    x(i+1) = x(i) + time_step*x_dot(i);
    fun = @(X2_dot)(d3*(x_dot(i) - X2_dot)^3 + d1*(x_dot(i) - X2_dot) - k3*x2(i)^3 - k1*x2(i));
%     fun = @(X2_dot)(d3*(x_dot(i+1) - X2_dot)^3 + d1*(x_dot(i+1) - X2_dot) - k3*x2(i+1)^3 - k1*x2(i+1));

    x2_dot(i+1) = fzero(fun, x2_dot(i));
    
    x2(i+1) = x2(i) + time_step*x2_dot(i);
    F(i+1)=k3*x2(i+1)^3+k1*x2(i+1);
end

% Plot Force v Extension

% extension = x2;
% force = k3*(x2-x1).^3 + k1*(x2-x1);
%force = d*v1;

Fs = k3*x2.^3+k1*x2;
Fd = d3*(x_dot - x2_dot).^3 + d1*(x_dot - x2_dot);
F = 0.5*(Fs+Fd);
figure; plot(Fs); hold on; plot(Fd); hold off;

% figure; plot(extension,force); %hold on; plot(t,force); legend('F-L','F(t)');
% figure; plot(x(20/time_step:end),F(20/time_step:end));
% figure; plot(x,F);
