% Alternative Kelvin Model

% Kelvin model

time_step = 0.01;
end_time = 100;
t = linspace(0,end_time,end_time/time_step + 1);
x = zeros(size(t));
x_dot = zeros(size(t));
F = zeros(size(t));
x2 = zeros(size(t));
x2_dot = zeros(size(t));
x1_dot = zeros(size(t));


p=1; B=1;
ks = 20; kp=2; d=10;
ks1 = 20; ks3 = 10;
kp1 = 5; kp3 = 7;
d1 = 10; d3 =4;

% sinusoidal forcing
f_t = B*sin(p*t);
% step forcing
% f_t = ones(size(t));


% Euler Integration

for i=1:length(t)-1
%     F(i+1) = F(i) + time_step*(ks/d*(kp*x(i)+d*x_dot(i)-(1+kp/ks)*F(i) ));
    
    x_dot(i+1) = x_dot(i) + time_step*(f_t(i)-F(i));
    x(i+1) = x(i) + time_step*x_dot(i);
    
    fun = @(X2_dot)(d3*(x_dot(i) - X2_dot)^3 + d1*(x_dot(i) - X2_dot) + kp3*(x(i) - x2(i))^3 + kp1*(x(i) - x2(i)) - ks3*x2(i)^3 - ks1*x2(i));
    x2_dot(i+1) = fzero(fun, x2_dot(i));

    x2(i+1) = x2(i) + time_step*x2_dot(i);
    x1_dot(i+1) = x_dot(i+1)-x2_dot(i+1);
    
    Fsee = ks1*x2(i+1)+ks3*x2(i+1)^3;
    Fsde = d1*(x_dot(i+1) - x2_dot(i+1)) + d3*(x_dot(i+1) - x2_dot(i+1))^3;
    Fpee = kp1*x(i+1) + kp3*x(i+1)^3;
    
    
    
    F(i+1) = 0.5*(Fsee+Fsde+Fpee);
    
    
end

% Plot Force v Extension

figure; plot(x(30/time_step:end),F(30/time_step:end)); xlabel('Extension (m)'); ylabel('Force (N)'); title('Nonlinear Kelvin Model B');
figure; plot(x,F); xlabel('Extension (m)'); ylabel('Force (N)'); title('Noninear Kelvin Model B');


