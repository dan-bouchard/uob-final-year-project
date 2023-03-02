% Voight Model

time_step = 0.01;
end_time = 100;
t = linspace(0,end_time,end_time/time_step + 1);

x = zeros(size(t));
x_dot = zeros(size(t));
F = zeros(size(t));

p=2; B=1; k=5; d=10;
k1 = 20; k3 = 10;
d1 = 10; d3 = 20;


% sinusoidal forcing
f_t = B*sin(p*t);
% step forcing
% f_t = ones(size(t));

% Euler Integration

for i=1:length(t)-1    
    x_dot(i+1) = x_dot(i) + time_step*(f_t(i)-F(i));
    x(i+1) = x(i) + time_step*x_dot(i);
    
    F(i+1)=k1*x(i+1) + k3*x(i+1)^3 + d1*x_dot(i+1) + d3*x_dot(i+1)^3;
%     F(i+1)=k1*x(i+1);
end

figure; plot(x(20/time_step:end),F(20/time_step:end)); xlabel('Extension (m)'); ylabel('Force (N)'); title('Nonlinear Voight Model');
figure; plot(t,x); 